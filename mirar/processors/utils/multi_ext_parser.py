"""
Module for loading and creating multi-extension FITS (MEF) images
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
from mirar.errors import ImageNotFoundError
from mirar.io import open_fits
from mirar.paths import (
    RAW_IMG_SUB_DIR,
    base_output_dir,
    base_raw_dir,
    get_output_dir,
    get_output_path,
)
from mirar.processors.base_processor import BaseImageProcessor
from mirar.processors.utils.image_loader import unzip

logger = logging.getLogger(__name__)


class MultiExtParser(BaseImageProcessor):
    """Processor to split multi-extension FITS (MEF) images into single-extension files.
    Should be run before ImageLoader, especially for SEDMv2
    """

    base_key = "load"  # should this be changed?

    def __init__(  # pylint: disable=too-many-arguments
        self,
        input_sub_dir: str = RAW_IMG_SUB_DIR,
        output_sub_dir: str = "raw_split",
        input_img_dir: str = base_raw_dir,
        output_img_dir: str | Path = base_output_dir,
        load_image: Callable[[str], [np.ndarray, astropy.io.fits.Header]] = open_fits,
        skip_first: bool = False,
        extension_num_header_key: str = None,
        only_extract_num: int = None,
    ):
        """
        :param input_sub_dir: subdirectory to look for images
        :param output_sub_dir: subdirectory to save split single extension images
        :param input_img_dir: parent directory of input_sub_dir
        :param output_sub_dir: parent directory of output_sub_dir
        :param load_image: function to load image
        :param skip_first: option to skip the first frame in MEF
        :param extension_num_header_key: If provided, will use the corresponding value
        in the header to identify an extension_number for every image and save the file
        as <>_{extension_number}.fits. If None, will serially number the extensions.
        :param only_extract_num: If provided, will only extract the extension with this
        number. If None, will extract all extensions. extension_number is calculated
        as described in extension_num_header_key.
        """
        super().__init__()
        self.input_sub_dir = input_sub_dir
        self.input_img_dir = input_img_dir
        self.output_img_dir = output_img_dir
        self.load_image = load_image
        self.skip_first = skip_first
        self.output_sub_dir = output_sub_dir
        self.extension_num_header_key = extension_num_header_key
        self.only_extract_num = only_extract_num

        if isinstance(self.output_img_dir, str):
            self.output_img_dir = Path(self.output_img_dir)

    def __str__(self):
        return (
            f"Processor to parse MEF images from the {self.input_sub_dir} subdirectory "
        )

    def parse(self, path: str) -> list:
        """
        Function to open a raw MEF image, write each extension to a new file

        :param path: path of raw MEF image
        :return: new paths of single-extension files
        """

        output_dir = get_output_dir(
            dir_root=self.output_sub_dir,
            sub_dir=self.night_sub_dir,
            output_dir=self.output_img_dir,
        )
        if not output_dir.exists():
            output_dir.mkdir(parents=True)

        new_paths = []
        with astropy.io.fits.open(path) as hdu:
            logger.debug(f"This file - {path} - has {len(hdu)} extensions.")

            hdr0 = hdu[0].header  # pylint: disable=no-member
            # zip hdr0's values and comments
            zipped = list(zip(hdr0.values(), hdr0.comments))
            # combining main header (hdr0) with extension header
            if self.skip_first:
                start = 2
                logger.debug("Ignoring first extension frame")
            else:
                start = 1
            for ext in range(start, len(hdu)):
                hdrext = hdu[ext].header

                extension_num_str = str(ext)
                if self.extension_num_header_key is not None:
                    extension_num_str = hdrext[self.extension_num_header_key]

                if self.only_extract_num is not None:
                    if int(extension_num_str) != self.only_extract_num:
                        continue

                # append hdr0 to hdrext
                for count, key in enumerate(list(hdr0.keys())):
                    hdrext.append((key, zipped[count][0], zipped[count][1]))

                # save to new file with 1 extension
                splitfile_basename = (
                    f"{os.path.basename(path).split('.fits')[0]}_"
                    f"{extension_num_str}.fits"
                )

                splitfile_path = get_output_path(
                    base_name=splitfile_basename,
                    dir_root=self.output_sub_dir,
                    sub_dir=self.night_sub_dir,
                    output_dir=self.output_img_dir,
                )

                try:
                    astropy.io.fits.writeto(
                        splitfile_path, hdu[ext].data, hdrext, overwrite=True
                    )  # pylint: disable=no-member
                    new_paths.append(splitfile_path)
                except TypeError:
                    logger.error(
                        f"Failed to write {splitfile_path} to file. "
                        f"This image is likely corrupted."
                    )

        return new_paths

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        input_dir = os.path.join(
            self.input_img_dir, os.path.join(self.night_sub_dir, self.input_sub_dir)
        )
        return load_from_dir(input_dir, parse_f=self.parse)


def load_from_dir(
    input_dir: str | Path, parse_f: Callable[[list[str | Path]], Image]
) -> ImageBatch:
    """
    Function to parse all MEF images in a directory

    :param input_dir: directory path
    :param parse_f: function to parse MEF images
    :return: nothing...
    """

    img_list = sorted(glob(f"{input_dir}/*.fits"))

    # check for zipped files too
    zipped_list = sorted(glob(f"{input_dir}/*.fz"))
    if len(zipped_list) > 0:
        unzipped_list = unzip(zipped_list)
        for file in unzipped_list:
            img_list.append(file)

    logger.info(f"Loading from {input_dir}, with {len(img_list)} images")

    if len(img_list) < 1:
        err = f"No images found in {input_dir}. Please check path is correct!"
        logger.error(err)
        raise ImageNotFoundError(err)

    for path in tqdm(img_list):
        parse_f(path)

    empty_batch = ImageBatch()
    return empty_batch
