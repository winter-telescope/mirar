import copy
import logging

import astropy.io.fits
import numpy as np

from winterdrp.data import Dataset, Image, ImageBatch
from winterdrp.paths import BASE_NAME_KEY
from winterdrp.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

sub_id_key = "SUBID"


class SplitImage(BaseImageProcessor):

    base_key = "split"

    def __init__(
        self, buffer_pixels: int = 0, n_x: int = 1, n_y: int = 1, *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.buffer_pixels = buffer_pixels
        self.n_x = n_x
        self.n_y = n_y

    def __str__(self) -> str:
        return f"Processor to split images into {self.n_x*self.n_y} smaller images."

    def get_range(
        self,
        n: int,
        pixel_width: int,
        i: int,
    ) -> tuple[int, int]:
        lower = max(0, i * int(pixel_width / n) - self.buffer_pixels)
        upper = min(pixel_width, (1 + i) * int(pixel_width / n) + self.buffer_pixels)
        return lower, upper

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        new_images = ImageBatch()

        logger.info(f"Splitting each data into {self.n_x*self.n_y} sub-images")

        for i, image in enumerate(batch):

            pix_width_x, pix_width_y = image.get_data().shape

            k = 0

            for ix in range(self.n_x):

                x_0, x_1 = self.get_range(self.n_x, pix_width_x, ix)

                for iy in range(self.n_y):
                    y_0, y_1 = self.get_range(self.n_y, pix_width_y, iy)

                    new_data = np.array(image.get_data()[x_0:x_1, y_0:y_1])

                    new_header = copy.copy(image.get_header())

                    for key in ["DETSIZE", "INFOSEC", "TRIMSEC", "DATASEC"]:
                        if key in new_header.keys():
                            del new_header[key]

                    sub_img_id = f"{ix}_{iy}"

                    new_header["SUBCOORD"] = (
                        sub_img_id,
                        "Sub-data coordinate, in form x_y",
                    )

                    new_header[sub_id_key] = k
                    k += 1

                    new_header["SRCIMAGE"] = (
                        image[BASE_NAME_KEY],
                        "Source data name, from which sub-data was made",
                    )

                    new_header["NAXIS1"], new_header["NAXIS2"] = new_data.shape

                    new_header[BASE_NAME_KEY] = image[BASE_NAME_KEY].replace(
                        ".fits", f"_{sub_img_id}.fits"
                    )
                    new_images.append(Image(data=new_data, header=new_header))

        return new_images

    def update_dataset(self, dataset: Dataset) -> Dataset:

        all_new_batches = []

        for batch in dataset:
            new_batches = [[] for _ in range(self.n_x * self.n_y)]

            for i, image in enumerate(batch):
                idx = image[sub_id_key]
                new_batches[idx] += [image]

            all_new_batches += new_batches

        return Dataset(all_new_batches)
