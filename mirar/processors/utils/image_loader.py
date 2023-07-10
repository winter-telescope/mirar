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
from mirar.errors import ImageNotFoundError, ProcessorError
from mirar.io import (
    MissingCoreFieldError,
    check_file_is_complete,
    check_image_has_core_fields,
    open_fits,
)
from mirar.paths import RAW_IMG_KEY, RAW_IMG_SUB_DIR, base_raw_dir
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class BadImageError(ProcessorError):
    """Exception for bad images"""


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
            f"Processor to load images from the '{self.input_sub_dir}' subdirectory "
            f"using the '{self.load_image.__name__}' function"
        )

    def open_raw_image(self, path: str) -> Image:
        """
        Function to open a raw image as an Image object

        :param path: path of raw image
        :return: Image object
        """
        data, header = self.load_image(path)

        new_img = Image(data.astype(np.float64), header)

        try:
            check_image_has_core_fields(new_img)
        except MissingCoreFieldError as err:
            raise BadImageError from err

        return new_img

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

    :param input_dir: Input directory
    :param open_f: Function to open images
    :return: ImageBatch object
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
        if check_file_is_complete(path):
            image = open_f(path)
            images.append(image)
        else:
            logger.warning(f"File {path} is not complete. Skipping!")

    return images


def unzip(zipped_list: list[str]) -> list[str]:
    """
    Function to unzip a list of files?

    :param zipped_list: List of zipped files
    :return: List of renamed files
    """
    unzipped_list = [file.split(".fz")[0] for file in zipped_list]
    for i, file in enumerate(zipped_list):
        os.rename(file, unzipped_list[i])

    return unzipped_list


class LoadImageFromHeader(BaseImageProcessor):
    """
    Class to load images from header information
    """

    base_key = "load_from_header"

    def __init__(
        self,
        header_key: str = RAW_IMG_KEY,
        copy_header_keys: str | list[str] = None,
        load_image: Callable[[str], [np.ndarray, astropy.io.fits.Header]] = open_fits,
    ):
        super().__init__()
        self.header_key = header_key
        self.copy_header_keys = copy_header_keys
        self.load_image = load_image
        if isinstance(self.copy_header_keys, str):
            self.copy_header_keys = [self.copy_header_keys]

    def __str__(self):
        return f"Processor to load images from header key {self.header_key}"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        new_batch = ImageBatch()
        for image in batch:
            new_image_file = image.header[self.header_key]
            new_image_data, new_header = self.load_image(new_image_file)
            new_image = Image(new_image_data, new_header)
            if self.copy_header_keys is not None:
                for key in self.copy_header_keys:
                    new_image.header[key] = image.header[key]

            try:
                check_image_has_core_fields(new_image)
            except MissingCoreFieldError as err:
                raise BadImageError from err

            new_batch.append(new_image)

        return new_batch
