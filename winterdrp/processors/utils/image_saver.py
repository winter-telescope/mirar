import os

import astropy.io.fits
import numpy as np
from winterdrp.paths import output_path, output_dir
from winterdrp.processors.base_processor import BaseProcessor

latest_save_key = "SAVEPATH"


class ImageSaver(BaseProcessor):

    base_key = "save"

    def __init__(
            self,
            output_dir_name: str,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.output_dir_name = output_dir_name

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        try:
            os.makedirs(output_dir(dir_root=self.output_dir_name, sub_dir=self.night_sub_dir))
        except OSError:
            pass

        for i, img in enumerate(images):

            header = headers[i]

            path = output_path(
                header["BASENAME"],
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir
            )

            self.save_fits(img, header, path)
            header[latest_save_key] = path

        return images, headers
