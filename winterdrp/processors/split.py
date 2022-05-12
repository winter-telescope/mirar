import copy

import astropy.io.fits
import numpy as np
import logging
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.paths import base_name_key

logger = logging.getLogger(__name__)


class SplitImage(BaseProcessor):

    base_key = "split"

    def __init__(
            self,
            buffer_pixels: int = 0,
            n_x: int = 1,
            n_y: int = 1,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.buffer_pixels = buffer_pixels
        self.n_x = n_x
        self.n_y = n_y

    def get_range(
            self,
            n: int,
            pixel_width: int,
            i: int,
    ) -> tuple[int, int]:
        lower = max(0, i * int(pixel_width/n) - self.buffer_pixels)
        upper = min(pixel_width, (1+i)*int(pixel_width/n) + self.buffer_pixels)
        return lower, upper

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        new_images = []
        new_headers = []

        logger.info(f"Splitting each image into {self.n_x*self.n_y} sub-images")

        for i, data in enumerate(images):

            base_header = headers[i]

            pix_width_x, pix_width_y = data.shape

            k = 0

            for ix in range(self.n_x):

                x_0, x_1 = self.get_range(self.n_x, pix_width_x, ix)

                for iy in range(self.n_y):
                    y_0, y_1 = self.get_range(self.n_y, pix_width_y, iy)

                    new_data = np.array(data[x_0:x_1, y_0:y_1])

                    new_header = copy.copy(base_header)

                    for key in ["DETSIZE", "INFOSEC", "TRIMSEC", "DATASEC"]:
                        if key in new_header.keys():
                            del new_header[key]

                    sub_img_id = f"{ix}_{iy}"

                    new_header["SUBCOORD"] = (sub_img_id, "Sub-image coordinate, in form x_y")

                    new_header["SUBID"] = k
                    k += 1

                    new_header["SRCIMAGE"] = (
                        base_header[base_name_key],
                        "Source image name, from which sub-image was made"
                    )

                    new_header[base_name_key] = base_header[base_name_key].replace(
                        ".fits", f"_{sub_img_id}.fits"
                    )
                    new_images.append(new_data)
                    new_headers.append(new_header)

        return new_images, new_headers
