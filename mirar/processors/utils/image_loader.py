"""
Module for loading images
"""

import logging
import os
from collections.abc import Callable
from glob import glob
from pathlib import Path

from tqdm import tqdm

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError, NoncriticalProcessingError, ProcessorError
from mirar.io import (
    MissingCoreFieldError,
    check_file_is_complete,
    check_image_has_core_fields,
    open_mef_image,
    open_raw_image,
)
from mirar.paths import RAW_IMG_KEY, RAW_IMG_SUB_DIR, base_raw_dir
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class BadImageError(ProcessorError):
    """Exception for bad images"""


class InvalidImage(NoncriticalProcessingError):
    """Image should be skipped"""


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


def load_from_list(
    img_list: list[str | Path],
    open_f: Callable[[str | Path], Image | list[Image]],
) -> ImageBatch:
    """
    Load images from a list of files

    :param img_list: Image list
    :param open_f: Function to open images
    :return: ImageBatch object
    """
    images = ImageBatch()

    for path in tqdm(img_list):
        if check_file_is_complete(path):
            try:
                image_list = open_f(path)

                if not isinstance(image_list, list):
                    image_list = [image_list]

                for image in image_list:
                    try:
                        check_image_has_core_fields(image)
                    except MissingCoreFieldError as err:
                        raise BadImageError(err) from err
                    images.append(image)
            except InvalidImage:
                logger.warning(f"Image {path} is invalid. Skipping!")
            except BadImageError:
                logger.error(f"Image {path} cannot be parsed. Skipping!")
        else:
            logger.warning(f"File {path} is not complete. Skipping!")

    return images


def load_from_dir(
    input_dir: str | Path,
    open_f: Callable[[str | Path], Image | list[Image]],
    bad_image_file: Path | None = None,
) -> ImageBatch:
    """
    Function to load all images in a directory

    :param input_dir: Input directory
    :param open_f: Function to open images
    :param bad_image_file: File containing list of bad images to skip
    :return: ImageBatch object
    """
    img_list = sorted(glob(f"{input_dir}/*.fits"))

    # check for zipped files too
    zipped_list = sorted(glob(f"{input_dir}/*.fz"))
    if len(zipped_list) > 0:
        unzipped_list = unzip(zipped_list)
        for file in unzipped_list:
            img_list.append(file)

    if bad_image_file is not None:

        with open(bad_image_file, "rb") as f:
            bad_images = f.read().decode().splitlines()

        clean_list = [x for x in img_list if os.path.basename(x) not in bad_images]

        n_skipped = len(img_list) - len(clean_list)

        logger.info(f"Skipping {n_skipped} bad images listed in {bad_image_file}")

        img_list = clean_list

    if len(img_list) < 1:
        err = f"No images found in {input_dir}. Please check path is correct!"
        logger.error(err)
        raise ImageNotFoundError(err)

    return load_from_list(img_list, open_f)


class ImageLoader(BaseImageProcessor):
    """Processor to load raw images."""

    base_key = "load"

    image_type = Image
    default_load_image = staticmethod(open_raw_image)

    def __init__(
        self,
        input_sub_dir: str = RAW_IMG_SUB_DIR,
        input_img_dir: str | Path = base_raw_dir,
        load_image: Callable[[str], Image | list[Image]] = None,
        bad_image_file: str = "bad_images.txt",
    ):
        super().__init__()
        self.input_sub_dir = input_sub_dir
        self.input_img_dir = Path(input_img_dir)
        if load_image is None:
            load_image = self.default_load_image
        self.load_image = load_image
        self.bad_image_file = bad_image_file

    def description(self):
        return (
            f"Processor to load images from the '{self.input_sub_dir}' subdirectory "
            f"using the '{self.load_image.__name__}' function"
        )

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        input_dir = self.input_img_dir.joinpath(
            os.path.join(self.night_sub_dir, self.input_sub_dir)
        )
        bad_image_path = input_dir / self.bad_image_file

        bad_image_file = bad_image_path if bad_image_path.is_file() else None

        return load_from_dir(
            input_dir,
            open_f=self.load_image,
            bad_image_file=bad_image_file,
        )


class LoadImageFromHeader(BaseImageProcessor):
    """
    Class to load images from header information
    """

    base_key = "load_from_header"

    def __init__(
        self,
        header_key: str = RAW_IMG_KEY,
        copy_header_keys: str | list[str] = None,
        load_image: Callable[[str | Path], Image] = open_raw_image,
    ):
        super().__init__()
        self.header_key = header_key
        self.copy_header_keys = copy_header_keys
        self.load_image = load_image
        if isinstance(self.copy_header_keys, str):
            self.copy_header_keys = [self.copy_header_keys]

    def description(self):
        return f"Processor to load images from header key {self.header_key}"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        new_batch = ImageBatch()
        for image in batch:
            new_image_file = image.header[self.header_key]
            new_image = self.load_image(new_image_file)

            if self.copy_header_keys is not None:
                for key in self.copy_header_keys:
                    new_image.header[key] = image.header[key]

            try:
                check_image_has_core_fields(new_image)
            except MissingCoreFieldError as err:
                raise BadImageError(err) from err

            new_batch.append(new_image)

        return new_batch


class ImageListLoader(BaseImageProcessor):
    """Processor to load raw images."""

    base_key = "loadlist"

    image_type = Image
    default_load_image = staticmethod(open_raw_image)

    def __init__(
        self,
        img_list: list[Path],
        load_image: Callable[[str], Image | list[Image]] = None,
    ):
        super().__init__()
        self.img_list = img_list
        if len(self.img_list) < 1:
            err = "No images found in list. Please check path is correct!"
            logger.error(err)
            raise ImageNotFoundError(err)
        if load_image is None:
            load_image = self.default_load_image
        self.load_image = load_image

    def description(self):
        return (
            f"Load {len(self.img_list)} images from "
            f"list of {len(self.img_list)} files"
        )

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        """
        Load images from a list of files

        :param batch: Batch of images
        :return: New batch of images
        """

        if len(batch) > 0:
            logger.warning("Batch is not empty. Overwriting images!")

        return load_from_list(
            self.img_list,
            open_f=self.load_image,
        )


class MEFLoader(ImageLoader):
    """Processor to load MEF images."""

    base_key = "load_mef"
    default_load_image = staticmethod(open_mef_image)
