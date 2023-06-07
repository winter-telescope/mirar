"""
Module with classes to perform photometry on an image or candidates
"""
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

from mirar.data import Image
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_SAVE_KEY,
    NORM_PSFEX_KEY,
    UNC_IMG_KEY,
    XPOS_KEY,
    YPOS_KEY,
    get_output_dir,
)
from mirar.processors.base_processor import (
    BaseDataframeProcessor,
    BaseImageProcessor,
    BaseProcessor,
    ImageHandler,
)
from mirar.processors.photometry.utils import (
    aper_photometry,
    get_rms_image,
    make_cutouts,
    make_psf_shifted_array,
    psf_photometry,
)

logger = logging.getLogger(__name__)


class BasePhotometry:
    """
    Parent class to run photometry on cutouts
    """

    def __init__(self):
        pass

    def perform_photometry(self, image_cutout: np.array, unc_image_cutout: np.array):
        """
        Function to perform photometry
        Args:
            image_cutout: 2D image cutout array
            unc_image_cutout: 2D image uncertainty cutout

        Returns:

        """
        raise NotImplementedError


class PSFPhotometry(BasePhotometry):
    """
    Class to run PSF photometry on cutouts
    """

    def __init__(self, psf_filename: str):
        super().__init__()
        self.psf_filename = psf_filename

    def perform_photometry(self, image_cutout: np.array, unc_image_cutout: np.array):
        psfmodels = make_psf_shifted_array(
            psf_filename=self.psf_filename,
            cutout_size_psf_phot=int(image_cutout.shape[0] / 2),
        )

        flux, fluxunc, minchi2, xshift, yshift = psf_photometry(
            image_cutout, unc_image_cutout, psfmodels
        )
        return flux, fluxunc, minchi2, xshift, yshift


class AperturePhotometry(BasePhotometry):
    """
    Class to run aperture photometry on cutouts
    """

    def __init__(
        self,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
    ):
        if not isinstance(aper_diameters, list):
            aper_diameters = [aper_diameters]
        if not isinstance(bkg_in_diameters, list):
            bkg_in_diameters = [bkg_in_diameters]
        if not isinstance(bkg_out_diameters, list):
            bkg_out_diameters = [bkg_out_diameters]
        super().__init__()

        self.aper_diameters = aper_diameters
        self.bkg_in_diameters = bkg_in_diameters
        self.bkg_out_diameters = bkg_out_diameters

    def perform_photometry(self, image_cutout: np.array, unc_image_cutout: np.array):
        fluxes, fluxuncs = [], []
        for ind, aper_diam in enumerate(self.aper_diameters):
            flux, fluxunc = aper_photometry(
                image_cutout,
                unc_image_cutout,
                aper_diam,
                self.bkg_in_diameters[ind],
                self.bkg_out_diameters[ind],
            )
            fluxes.append(flux)
            fluxuncs.append(fluxunc)
        return fluxes, fluxuncs


class BasePhotometryProcessor(BaseProcessor, ImageHandler):
    """
    Parent processor to run photometry
    """

    def __init__(
        self, phot_cutout_size: int = 20, temp_output_sub_dir: str = "photometry"
    ):
        super().__init__()
        self.phot_cutout_size = phot_cutout_size
        self.temp_output_sub_dir = temp_output_sub_dir
        self.photometry_out_temp_dir = get_output_dir(
            self.temp_output_sub_dir, self.night_sub_dir
        )

    def get_image_filenames(self, data_item: Image | pd.Series):
        """
        Get image and uncertainty image filenames
        Args:
            data_item: pandas DataFrame Series or astropy fits Header

        Returns:
            tuple: image name, uncertainty image name
        """
        raise NotImplementedError

    def get_physical_coordinates(self, data_item: Image | pd.Series):
        """
        Get physical coordinates on which to perform photometry
        Args:
            data_item: pandas DataFrame Series or astropy fits Header

        Returns:
            tuple: x, y coordinates at which photometry is performed

        """
        raise NotImplementedError

    def get_temp_image_unc_imagenames(self, image: Image) -> (str, str):
        self.photometry_out_temp_dir.mkdir(parents=True, exist_ok=True)

        image_basename = image.header[BASE_NAME_KEY]
        temp_imagepath = self.photometry_out_temp_dir.joinpath(image_basename)
        self.save_fits(image, temp_imagepath)

        unc_exists, temp_unc_imagepath = False, None
        if UNC_IMG_KEY in image.header.keys():
            unc_filename = image.header[UNC_IMG_KEY]
            if unc_filename is None:
                unc_exists = False
            else:
                temp_unc_imagepath = Path(unc_filename)
                unc_exists = os.path.exists(unc_filename)
        if not unc_exists:
            rms_image = get_rms_image(image)
            unc_filename = image[BASE_NAME_KEY] + ".unc"
            temp_unc_imagepath = self.photometry_out_temp_dir.joinpath(unc_filename)
            self.save_fits(rms_image, path=temp_unc_imagepath)
            logger.info(f"Saved unc file to {temp_unc_imagepath}")
        return temp_imagepath, temp_unc_imagepath

    def generate_cutouts(self, data_item: Image | pd.Series):
        """
        Generate image cutouts
        Args:
            data_item: pandas DataFrame Series or astropy fits Header

        Returns:
            tuple: 2D numpy arrays of the image cutout and uncertainty image cutout
        """
        imagename, unc_imagename = self.get_image_filenames(data_item)
        x, y = self.get_physical_coordinates(data_item)
        image_cutout, unc_image_cutout = make_cutouts(
            image_paths=[imagename, unc_imagename],
            position=(x, y),
            half_size=self.phot_cutout_size,
        )
        return image_cutout, unc_image_cutout


class BaseImagePhotometry(BasePhotometryProcessor, BaseImageProcessor):
    """
    Processor to run photometry on an image
    """

    def __init__(
        self,
        target_ra_key: str = "TARGRA",
        target_dec_key: str = "TARGDEC",
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.target_ra_key = target_ra_key
        self.target_dec_key = target_dec_key

    def get_image_filenames(self, image: Image):
        image_filename, unc_filename = self.get_temp_image_unc_imagenames(image)
        return image_filename, unc_filename

    def get_physical_coordinates(self, image: Image) -> (int, int):
        ra, dec = image[self.target_ra_key], image[self.target_dec_key]
        wcs = WCS(image.header)
        x, y = wcs.all_world2pix(ra, dec, 0)
        return int(np.round(x)), int(np.round(y))


class BaseCandidatePhotometry(BasePhotometryProcessor, BaseDataframeProcessor):
    """
    Processor to run photometry on a candidates table
    """

    def __init__(
        self,
        image_colname=LATEST_SAVE_KEY,
        unc_image_colname=UNC_IMG_KEY,
        psf_file_colname=NORM_PSFEX_KEY,
        x_colname=XPOS_KEY,
        y_colname=YPOS_KEY,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.image_key = image_colname
        self.unc_image_key = unc_image_colname
        self.psf_file_key = psf_file_colname
        self.xpos_key = x_colname
        self.ypos_key = y_colname

    def get_image_filenames(self, row):
        imagename = row[self.image_key]
        image = Image(header=fits.getheader(imagename), data=fits.getdata(imagename))

        image_filename, unc_filename = self.get_temp_image_unc_imagenames(image)

        return image_filename, unc_filename

    def get_physical_coordinates(self, row: pd.Series):
        x, y = row[self.xpos_key], row[self.ypos_key]
        return int(x), int(y)
