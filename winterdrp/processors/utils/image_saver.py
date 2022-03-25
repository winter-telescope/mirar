import os

import astropy.io.fits
import numpy as np
from winterdrp.paths import get_output_path, get_output_dir, latest_save_key, latest_mask_save_key
from winterdrp.processors.base_processor import BaseProcessor


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

            path = get_output_path(
                header["BASENAME"],
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir
            )

            header[latest_save_key] = path

            if self.write_mask:
                mask_path = self.save_mask(img, header, img_path=path)
                header[latest_mask_save_key] = mask_path

            self.save_fits(img, header, path)

        return images, headers
