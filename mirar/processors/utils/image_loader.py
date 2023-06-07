"""
Module for loading images
"""
import logging
import os
from collections.abc import Callable
from glob import glob
from pathlib import Path

import astropy.io.fits
import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.io import open_fits
from mirar.paths import RAW_IMG_SUB_DIR, base_raw_dir, core_fields
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class ImageLoader(BaseImageProcessor):
    """Processor to load raw images."""

    base_key = "load"

    def __init__(
        self,
        input_sub_dir: str = RAW_IMG_SUB_DIR,
        input_img_dir: str = base_raw_dir,
        load_image: Callable[[str], [np.ndarray, astropy.io.fits.Header]] = open_fits,
    ):
        super().__init__()
        self.input_sub_dir = input_sub_dir
        self.load_image = load_image
        self.input_img_dir = input_img_dir

    def __str__(self):
        return (
            f"Processor to load images from the {self.input_sub_dir} subdirectory "
            f"using the '{self.load_image.__name__}' function"
        )

    def open_raw_image(self, path: str) -> Image:
        """
        Function to open a raw image as an Image object

        :param path: path of raw image
        :return: Image object
        """
        data, header = self.load_image(path)

        for key in core_fields:
            if key not in header.keys():
                err = (
                    f"Essential key {key} not found in header. "
                    f"Please add this field first. Available fields are: "
                    f"{list(header.keys())}"
                )
                logger.error(err)
                raise KeyError(err)

        return Image(data.astype(np.float64), header)

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        input_dir = os.path.join(
            self.input_img_dir, os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        return load_from_dir(input_dir, open_f=self.open_raw_image)


def load_from_dir(
    input_dir: str | Path, open_f: Callable[[str | Path], Image]
) -> ImageBatch:
    """
    Function to load all images in a directory

    :param input_dir:
    :param open_f:
    :return:
    """

    img_list = sorted(glob(f"{input_dir}/*.fits"))

    # check for zipped files too
    zipped_list = sorted(glob(f"{input_dir}/*.fz"))
    if len(zipped_list) > 0:
        unzipped_list = unzip(zipped_list)
        for file in unzipped_list:
            img_list.append(file)

    if len(img_list) < 1:
        err = f"No images found in {input_dir}. Please check path is correct!"
        logger.error(err)
        raise ImageNotFoundError(err)

    images = ImageBatch()

    for path in img_list:
        image = open_f(path)
        images.append(image)

    return images


def unzip(zipped_list: list[str]) -> list[str]:
    unzipped_list = [file.split(".fz")[0] for file in zipped_list]
    for i, file in enumerate(zipped_list):
        os.rename(file, unzipped_list[i])

    return unzipped_list
