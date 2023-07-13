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
from tqdm import tqdm

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError, ProcessorError
from mirar.io import (
    MissingCoreFieldError,
    check_file_is_complete,
    check_image_has_core_fields,
    combine_mef_extension_file_headers,
    open_fits,
    open_mef_fits,
)
from mirar.paths import (
    BASE_NAME_KEY,
    RAW_IMG_KEY,
    RAW_IMG_SUB_DIR,
    base_raw_dir,
    core_fields,
    get_output_dir,
    get_temp_path,
)
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class BadImageError(ProcessorError):
    """Exception for bad images"""


class ImageLoader(BaseImageProcessor):
    """Processor to load raw images."""

    base_key = "load"

    image_type = Image

    def __init__(
        self,
        input_sub_dir: str = RAW_IMG_SUB_DIR,
        input_img_dir: str = base_raw_dir,
        load_image: Callable[[str], [np.ndarray, astropy.io.fits.Header]] = open_fits,
    ):
        super().__init__()
        self.input_sub_dir = input_sub_dir
        self.input_img_dir = input_img_dir
        self.load_image = load_image

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
            raise BadImageError(err) from err

        return new_img

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        input_dir = os.path.join(
            self.input_img_dir, os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        return load_from_dir(
            input_dir,
            open_f=self.open_raw_image,
        )


def load_from_dir(
    input_dir: str | Path,
    open_f: Callable[[str | Path], Image] | Callable[[str | Path], list[Image]],
    mef: bool = False,
) -> ImageBatch:
    """
    Function to load all images in a directory

    :param input_dir: Input directory
    :param open_f: Function to open images
    :param mef: Is this a mef frame?
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

    for path in tqdm(img_list):
        if check_file_is_complete(path):
            if not mef:
                image = open_f(path)
                images.append(image)
            else:
                image_list = open_f(path)
                for x in image_list:
                    images.append(x)
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
                raise BadImageError(err) from err

            new_batch.append(new_image)

        return new_batch


class MEFLoader(BaseImageProcessor):
    """Processor to load MEF images."""

    base_key = "load_mef"

    def __init__(
        self,
        input_sub_dir: str = RAW_IMG_SUB_DIR,
        input_img_dir: str = base_raw_dir,
        load_image: Callable[
            [str],
            [astropy.io.fits.Header, list[np.ndarray], list[astropy.io.fits.Header]],
        ] = open_mef_fits,
        extension_num_header_key: str = None,
        only_extract_num: int = None,
    ):
        super().__init__()
        self.input_sub_dir = input_sub_dir
        self.input_img_dir = input_img_dir
        self.load_image = load_image
        self.extension_num_header_key = extension_num_header_key
        self.only_extract_num = only_extract_num

    def open_split_raw_image(self, path: str | Path) -> list[Image]:
        """
        Function to open a raw image as an Image object

        :param path: path of raw image
        :return: Image object
        """
        primary_header, ext_data_list, ext_header_list = self.load_image(path)

        for key in core_fields:
            if key not in primary_header.keys():
                err = (
                    f"Essential key {key} not found in header. "
                    f"Please add this field first. Available fields are: "
                    f"{list(primary_header.keys())}"
                )
                logger.error(err)
                raise KeyError(err)

        ext_data_list = [x.astype(np.float64) for x in ext_data_list]
        split_images_list = []

        temp_dir = get_output_dir(dir_root="temp_load_mef", sub_dir=self.night_sub_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)

        temp_files = []
        for ext_num, ext_data in enumerate(ext_data_list):
            ext_header = ext_header_list[ext_num]

            extension_num_str = str(ext_num)

            if self.extension_num_header_key is not None:
                extension_num_str = ext_header[self.extension_num_header_key]

            if self.only_extract_num is not None:
                if int(extension_num_str) != self.only_extract_num:
                    continue

            # append primary_header to hdrext
            new_single_header = combine_mef_extension_file_headers(
                primary_header=primary_header, extension_header=ext_header
            )

            new_single_header[BASE_NAME_KEY] = (
                f"{primary_header[BASE_NAME_KEY].split('.fits')[0]}_"
                f"{extension_num_str}.fits"
            )

            # TODO : For some reason the above step doesn't gel well with astropy, and
            # it can't load the header. So we save it to a temp file and then load it,
            # which seems to work.
            temp_savepath = get_temp_path(
                output_dir=temp_dir, file_path=new_single_header[BASE_NAME_KEY]
            )
            astropy.io.fits.writeto(
                temp_savepath, ext_data, ext_header, overwrite=True
            )  # pylint: disable=no-member
            temp_files.append(temp_savepath)
            tmp_data, tmp_header = open_fits(temp_savepath)
            split_images_list.append(Image(data=tmp_data, header=tmp_header))

        for temp_file in temp_files:
            temp_file.unlink()

        return split_images_list

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        input_dir = os.path.join(
            self.input_img_dir, os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        return load_from_dir(input_dir, open_f=self.open_split_raw_image, mef=True)
