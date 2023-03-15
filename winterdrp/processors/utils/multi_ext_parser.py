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

from winterdrp.data import Image, ImageBatch
from winterdrp.errors import ImageNotFoundError
from winterdrp.paths import RAW_IMG_SUB_DIR, base_raw_dir
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.processors.utils import ImageLoader

logger = logging.getLogger(__name__)


class MultiExtParser(BaseImageProcessor):
    """Processor to split multi-extension FITS (MEF) images into single-extension files.
    Should be run before ImageLoader, especially for SEDMv2
    """

    base_key = "load"  # TODO: should this be changed?

    def __init__(
        self,
        input_sub_dir: str = RAW_IMG_SUB_DIR,
        input_img_dir: str = base_raw_dir,
    ):
        super().__init__()
        self.input_sub_dir = input_sub_dir
        self.input_img_dir = input_img_dir

    def __str__(self):
        return (
            f"Processor to parse MEF images from the {self.input_sub_dir} subdirectory "
        )

    def parse(self, path: str) -> Image:
        """
        Function to open a raw MEF image, write each extension to a new file

        :param path: path of raw MEF image
        :return: ????????/ ###Image objects (or ImageBatch?)

        ** should only run on science images **
        *** need to manually place MEF science images in a /mef/ subdirectory,
        ex: /[instrument]/[night]/raw/mef/ ***
        """

        with astropy.io.fits.open(path) as hdu:
            num_ext = len(hdu)
            print("This file has ", num_ext, "extensions")

            hdr0 = hdu[0].header
            # zip hdr0's values and comments
            zipped = list(zip(hdr0.values(), hdr0.comments))
            # combining main header (hdr0) with extension header
            for ext in range(1, num_ext):
                data = hdu[ext].data
                hdrext = hdu[ext].header

                # append hdr0 to hdrext
                for count, key in enumerate(list(hdr0.keys())):
                    hdrext.append((key, zipped[count][0], zipped[count][1]))

                # save to new file with 1 extension
                notmefpath = path.split("/mef/")[0] + path.split("/mef")[1]
                newpath = notmefpath.split(".fits")[0] + "_" + str(ext) + ".fits"
                astropy.io.fits.writeto(
                    newpath, data, hdrext, overwrite=True
                )  # pylint: disable=no-member

        msg = "parse function complete... now what?"
        return msg

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        input_dir = os.path.join(
            self.input_img_dir, os.path.join(self.night_sub_dir, self.input_sub_dir)
        )
        print("input_dir: ", input_dir)
        return load_from_dir(input_dir, open_f=self.parse)


def load_from_dir(
    input_dir: str | Path, open_f: Callable[[str | Path], Image]
) -> ImageBatch:
    """
    Function to load all images in a directory

    :param input_dir: directory path
    :param open_f: function to open raw images
    :return: ImageBatch
    """

    img_list = sorted(glob(f"{input_dir}/*.fits"))

    logger.info(f"Loading from {input_dir}, with {len(img_list)} images")

    if len(img_list) < 1:
        err = f"No images found in {input_dir}. Please check path is correct!"
        logger.error(err)
        raise ImageNotFoundError(err)

    images = ImageBatch()

    for path in img_list:
        image = open_f(path)
        images.append(image)

    return images