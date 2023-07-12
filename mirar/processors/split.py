"""
Module for splitting images into sub-images
"""
import copy
import logging

import numpy as np

from mirar.data import Dataset, Image, ImageBatch
from mirar.paths import BASE_NAME_KEY, LATEST_SAVE_KEY, LATEST_WEIGHT_SAVE_KEY
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

SUB_ID_KEY = "SUBID"
SUB_COORD_KEY = "SUBCOORD"


class SplitImage(BaseImageProcessor):
    """
    Processor for splitting images
    """

    base_key = "split"

    def __init__(self, buffer_pixels: int = 0, n_x: int = 1, n_y: int = 1):
        super().__init__()
        self.buffer_pixels = buffer_pixels
        self.n_x = n_x
        self.n_y = n_y

    def __str__(self) -> str:
        return f"Processor to split images into {self.n_x}x{self.n_y}={self.n_x*self.n_y} smaller images."

    def get_range(
        self,
        n_chunks: int,
        pixel_width: int,
        i: int,
    ) -> tuple[int, int]:
        """
        Function to return pixel index range for sub images

        :param n_chunks: number of chunks to divide axis into
        :param pixel_width: total pixel width of axis
        :param i: index of chunk to evaluate
        :return: lower pixel index and upper pixel index of chunk
        """
        lower = max(0, i * int(pixel_width / n_chunks) - self.buffer_pixels)
        upper = min(
            pixel_width, (1 + i) * int(pixel_width / n_chunks) + self.buffer_pixels
        )
        return lower, upper

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        new_images = ImageBatch()

        logger.debug(f"Splitting each data into {self.n_x*self.n_y} sub-images")

        for image in batch:
            pix_width_x, pix_width_y = image.get_data().shape

            k = 0

            for index_x in range(self.n_x):
                x_0, x_1 = self.get_range(self.n_x, pix_width_x, index_x)

                for index_y in range(self.n_y):
                    y_0, y_1 = self.get_range(self.n_y, pix_width_y, index_y)

                    new_data = np.array(image.get_data()[x_0:x_1, y_0:y_1])

                    new_header = copy.copy(image.get_header())

                    for key in ["DETSIZE", "INFOSEC", "TRIMSEC", "DATASEC"]:
                        if key in new_header.keys():
                            del new_header[key]

                    sub_img_id = f"{index_x}_{index_y}"

                    new_header[SUB_COORD_KEY] = (
                        sub_img_id,
                        "Sub-data coordinate, in form x_y",
                    )

                    new_header["SUBNX"] = (index_x + 1, "Sub-data x index")
                    new_header["SUBNY"] = (index_y + 1, "Sub-data y index")
                    new_header["SUBNXTOT"] = (self.n_x, "Total number of sub-data in x")
                    new_header["SUBNYTOT"] = (self.n_y, "Total number of sub-data in y")

                    new_header[SUB_ID_KEY] = k
                    k += 1

                    new_header["SRCIMAGE"] = (
                        image[BASE_NAME_KEY],
                        "Source data name, from which sub-data was made",
                    )

                    new_header["NAXIS1"], new_header["NAXIS2"] = new_data.shape

                    new_header[BASE_NAME_KEY] = image[BASE_NAME_KEY].replace(
                        ".fits", f"_{sub_img_id}.fits"
                    )

                    for key in [LATEST_SAVE_KEY, LATEST_WEIGHT_SAVE_KEY]:
                        if key in new_header.keys():
                            del new_header[key]

                    new_images.append(Image(data=new_data, header=new_header))

        return new_images

    def update_dataset(self, dataset: Dataset) -> Dataset:
        all_new_batches = []

        for batch in dataset:
            new_images = [[] for _ in range(self.n_x * self.n_y)]

            for image in batch:
                idx = image[SUB_ID_KEY]
                new_images[idx] += [image]

            all_new_batches += new_images
        all_new_batches = [ImageBatch(x) for x in all_new_batches]
        return Dataset(all_new_batches)
