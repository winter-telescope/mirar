"""
Module with classes to perform photometry on sources
"""

import logging
from abc import ABC
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from mirar.data import Image
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_SAVE_KEY,
    NORM_PSFEX_KEY,
    UNC_IMG_KEY,
    XPOS_KEY,
    YPOS_KEY,
    ZP_KEY,
    ZP_STD_KEY,
    get_output_dir,
)
from mirar.processors.base_processor import BaseSourceProcessor, ImageHandler
from mirar.processors.photometry.utils import get_rms_image, make_cutouts

logger = logging.getLogger(__name__)


class BasePhotometryProcessor(BaseSourceProcessor, ABC, ImageHandler):
    """
    Parent processor to run photometry
    """

    def __init__(
        self,
        phot_cutout_half_size: int = 20,
        temp_output_sub_dir: str = "photometry",
        zp_key: str = ZP_KEY,
        zp_std_key: str = ZP_STD_KEY,
        img_key: str = LATEST_SAVE_KEY,
        unc_image_key: str = UNC_IMG_KEY,
        psf_file_key: str = NORM_PSFEX_KEY,
        x_colname: str = XPOS_KEY,
        y_colname: str = YPOS_KEY,
        save_cutouts: bool = False,
    ):
        super().__init__()
        self.phot_cutout_half_size = phot_cutout_half_size
        self.temp_output_sub_dir = temp_output_sub_dir
        self.zp_key = zp_key
        self.zp_std_key = zp_std_key
        self.image_key = img_key
        self.unc_image_key = unc_image_key
        self.psf_file_key = psf_file_key
        self.xpos_key = x_colname
        self.ypos_key = y_colname
        self.save_cutouts = save_cutouts

    def save_temp_image(self, image) -> Path:
        """
        Save a temporary image and return its path

        :param image: Image object
        :return: Path to the temporary image
        """
        photometry_out_temp_dir = get_output_dir(
            self.temp_output_sub_dir, self.night_sub_dir
        )
        image_basename = image.header[BASE_NAME_KEY]
        temp_imagepath = photometry_out_temp_dir.joinpath(image_basename)
        temp_imagepath.parent.mkdir(parents=True, exist_ok=True)
        self.save_fits(image, temp_imagepath)
        return temp_imagepath

    def save_uncertainty_image(self, image: Image) -> Path:
        """
        Create an uncertainty image from the image and return the filenames of the
        image and uncertainty image

        :param image: Image object
        :return: Path to the uncertainty image
        """
        photometry_out_temp_dir = get_output_dir(
            self.temp_output_sub_dir, self.night_sub_dir
        )
        unc_filename = Path(image[BASE_NAME_KEY].replace(".fits", ".unc.fits"))
        rms_image = get_rms_image(image)
        unc_filename = photometry_out_temp_dir.joinpath(unc_filename)
        unc_filename.parent.mkdir(parents=True, exist_ok=True)
        self.save_fits(rms_image, path=unc_filename)
        logger.debug(f"Saved unc file to {unc_filename}")

        return unc_filename

    def generate_cutouts(
        self, imagename: Path, unc_imagename: Path, data_item: Image | pd.Series
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Generate image and uncertainty image cutouts. This function first saves the
        image and uncertainty image to temporary files, then generates the cutouts,
        then deletes the temporary files.
        :param imagename: Path to the image
        :param unc_imagename: Path to the uncertainty image
        :param data_item: pandas DataFrame Series or astropy fits Header

        :returns tuple: 2D numpy arrays of the image cutout and uncertainty image cutout
        """

        x, y = self.get_physical_coordinates(data_item)
        image_cutout, unc_image_cutout = make_cutouts(
            image_paths=[imagename, unc_imagename],
            position=(x, y),
            half_size=self.phot_cutout_half_size,
        )

        return image_cutout, unc_image_cutout

    def save_temp_image_uncimage(self, metadata: dict) -> tuple[Path, Path]:
        """
        Function to save the image and uncertainty image to temporary files

        :param metadata: Metadata dictionary
        :return: Tuple of image and uncertainty image filenames
        """
        imagename = metadata[self.image_key]

        image = Image(header=fits.getheader(imagename), data=fits.getdata(imagename))

        image_filename = self.save_temp_image(image)
        unc_filename = self.save_uncertainty_image(image)

        return image_filename, unc_filename

    def get_physical_coordinates(self, data_item: pd.Series) -> tuple[int, int]:
        """
        Get the physical coordinates of the source from the data item

        :param data_item: Series from the data table
        :return: X and Y coordinates of the source
        """
        row = data_item
        x, y = row[self.xpos_key], row[self.ypos_key]
        return int(x), int(y)
