"""
Module with classes to perform photometry on an image or candidates
"""
import logging
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
    ZP_KEY,
    ZP_STD_KEY,
    get_output_dir,
)
from mirar.processors.base_processor import (
    BaseImageProcessor,
    BaseProcessor,
    BaseSourceProcessor,
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

    def perform_photometry(
        self, image_cutout: np.array, unc_image_cutout: np.array
    ) -> tuple[float, float, float, float, float]:
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

    def perform_photometry(
        self, image_cutout: np.array, unc_image_cutout: np.array
    ) -> tuple[list[float], list[float]]:
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
        self,
        phot_cutout_size: int = 20,
        temp_output_sub_dir: str = "photometry",
        zp_key: str = ZP_KEY,
        zp_std_key: str = ZP_STD_KEY,
    ):
        """
        Args:
            phot_cutout_size: size of cutout to perform photometry on
            temp_output_sub_dir: subdirectory to save temporary photometry files
        """
        super().__init__()
        self.phot_cutout_size = phot_cutout_size
        self.temp_output_sub_dir = temp_output_sub_dir
        self.zp_key = zp_key
        self.zp_std_key = zp_std_key

    def save_temp_image_uncimage(
        self, data_item: Image | pd.Series
    ) -> tuple[Path, Path]:
        """
        Create temporary image and uncertainty images and return their paths.
        Args:
            data_item: pandas DataFrame Series or astropy fits Header

        Returns:
            tuple: image name, uncertainty image name
        """
        raise NotImplementedError

    def get_physical_coordinates(
        self, data_item: Image | pd.Series
    ) -> tuple[float, float]:
        """
        Get physical coordinates (x, y) at which to perform photometry
        Args:
            data_item: pandas DataFrame Series or astropy fits Header

        Returns:
            tuple: x, y coordinates at which photometry is performed

        """
        raise NotImplementedError

    def save_temp_image(self, image) -> Path:
        """
        Save a temporary image and return its path
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
        """
        photometry_out_temp_dir = get_output_dir(
            self.temp_output_sub_dir, self.night_sub_dir
        )
        unc_filename = Path(image[BASE_NAME_KEY] + ".unc")
        rms_image = get_rms_image(image)
        unc_filename = photometry_out_temp_dir.joinpath(unc_filename)
        unc_filename.parent.mkdir(parents=True, exist_ok=True)
        self.save_fits(rms_image, path=unc_filename)
        logger.debug(f"Saved unc file to {unc_filename}")

        return unc_filename

    def generate_cutouts(
        self, data_item: Image | pd.Series
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Generate image and uncertainty image cutouts. This function first saves the
        image and uncertainty image to temporary files, then generates the cutouts,
        then deletes the temporary files.
        Args:
            data_item: pandas DataFrame Series or astropy fits Header

        Returns:
            tuple: 2D numpy arrays of the image cutout and uncertainty image cutout
        """
        temp_imagename, temp_unc_imagename = self.save_temp_image_uncimage(data_item)
        x, y = self.get_physical_coordinates(data_item)
        image_cutout, unc_image_cutout = make_cutouts(
            image_paths=[temp_imagename, temp_unc_imagename],
            position=(x, y),
            half_size=self.phot_cutout_size,
        )
        temp_imagename.unlink()
        temp_unc_imagename.unlink()
        return image_cutout, unc_image_cutout


class BaseImagePhotometry(BasePhotometryProcessor, BaseImageProcessor):
    """
    Processor to run photometry on an image
    """

    def __init__(
        self,
        *args,
        target_ra_key: str = "TARGRA",
        target_dec_key: str = "TARGDEC",
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.target_ra_key = target_ra_key
        self.target_dec_key = target_dec_key

    def save_temp_image_uncimage(self, data_item: Image):
        image = data_item
        image_filename = self.save_temp_image(image)
        unc_filename = self.save_uncertainty_image(image)
        return image_filename, unc_filename

    def get_physical_coordinates(self, data_item: Image) -> (int, int):
        image = data_item
        ra, dec = image[self.target_ra_key], image[self.target_dec_key]
        wcs = WCS(image.header)
        x, y = wcs.all_world2pix(ra, dec, 0)
        return int(np.round(x)), int(np.round(y))


class BaseSourcePhotometry(BasePhotometryProcessor, BaseSourceProcessor):
    """
    Processor to run photometry on a candidates table
    """

    base_key = "base_source_photometry"

    def __init__(
        self,
        *args,
        image_colname=LATEST_SAVE_KEY,
        unc_image_colname=UNC_IMG_KEY,
        psf_file_colname=NORM_PSFEX_KEY,
        x_colname=XPOS_KEY,
        y_colname=YPOS_KEY,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.image_key = image_colname
        self.unc_image_key = unc_image_colname
        self.psf_file_key = psf_file_colname
        self.xpos_key = x_colname
        self.ypos_key = y_colname

    def save_temp_image_uncimage(self, data_item: pd.Series):
        row = data_item
        imagename = row[self.image_key]
        image = Image(header=fits.getheader(imagename), data=fits.getdata(imagename))

        image_filename = self.save_temp_image(image)
        unc_filename = self.save_uncertainty_image(image)

        return image_filename, unc_filename

    def get_physical_coordinates(self, data_item: pd.Series):
        row = data_item
        x, y = row[self.xpos_key], row[self.ypos_key]
        return int(x), int(y)
