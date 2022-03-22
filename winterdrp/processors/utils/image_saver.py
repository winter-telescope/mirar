import os

import astropy.io.fits
import numpy as np
from winterdrp.paths import output_path, get_output_dir, get_mask_path
from winterdrp.processors.base_processor import BaseProcessor

latest_save_key = "SAVEPATH"


class ImageSaver(BaseProcessor):

    base_key = "save"

    def __init__(
            self,
            output_dir_name: str,
            write_mask: bool = True,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.output_dir_name = output_dir_name
        self.write_mask = write_mask

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        try:
            os.makedirs(get_output_dir(dir_root=self.output_dir_name, sub_dir=self.night_sub_dir))
        except OSError:
            pass

        for i, img in enumerate(images):

            header = headers[i]

            path = output_path(
                header["BASENAME"],
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir
            )

            if self.write_mask:
                mask = np.isnan(img).astype(int)
                mask_path = get_mask_path(path)
                header["MASKPATH"] = mask_path
                self.save_fits(mask, header, mask_path)

            self.save_fits(img, header, path)
            header[latest_save_key] = path

        return images, headers
