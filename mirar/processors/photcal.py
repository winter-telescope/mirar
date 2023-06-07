"""
Module for running photometric calibration
"""
import logging
import os
from collections.abc import Callable
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats

from mirar.catalog.base_catalog import BaseCatalog
from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import copy_temp_file, get_output_dir, get_output_path
from mirar.processors.astromatic.sextractor.sextractor import (
    SEXTRACTOR_HEADER_KEY,
    Sextractor,
    sextractor_checkimg_map,
)
from mirar.processors.base_processor import BaseImageProcessor, PrerequisiteError
from mirar.processors.candidates.utils.regions_writer import write_regions_file
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)

# All the Sextractor parameters required for this script to run
REQUIRED_PARAMETERS = [
    "X_IMAGE",
    "Y_IMAGE",
    "FWHM_WORLD",
    "FLAGS",
    "ALPHAWIN_J2000",
    "DELTAWIN_J2000",
    "MAG_APER",
    "MAG_AUTO",
]


class PhotometryError(ProcessorError):
    """Base error for photometric calibration"""


class PhotometryReferenceError(PhotometryError):
    """Error related to the photometric reference catalogue"""


class PhotometrySourceError(PhotometryError):
    """Error related to the photometric source catalogue"""


class PhotometryCrossMatchError(PhotometryError):
    """Error related to cross-matching photometric reference and source catalogues"""


class PhotometryCalculationError(PhotometryError):
    """Error related to the photometric calibration"""


class PhotCalibrator(BaseImageProcessor):
    """
    Photometric calibrator processor
    """

    base_key = "photcalibrator"

    def __init__(
        self,
        ref_catalog_generator: Callable[[Image], BaseCatalog],
        temp_output_sub_dir: str = "phot",
        redo: bool = True,
        x_lower_limit: float = 100,
        x_upper_limit: float = 2800,  # Are these floats or ints?
        y_lower_limit: float = 100,
        y_upper_limit: float = 2800,
        fwhm_threshold_arcsec: float = 4.0,
        num_matches_threshold: int = 5,
        write_regions: bool = False,
        cache: bool = False,
    ):
        super().__init__()
        self.redo = redo  # What is this for?
        self.ref_catalog_generator = ref_catalog_generator
        self.temp_output_sub_dir = temp_output_sub_dir
        self.x_lower_limit = x_lower_limit
        self.x_upper_limit = x_upper_limit
        self.y_lower_limit = y_lower_limit
        self.y_upper_limit = y_upper_limit
        self.cache = cache

        # Why is this here not in catalog?
        self.fwhm_threshold_arcsec = fwhm_threshold_arcsec

        self.num_matches_threshold = num_matches_threshold
        self.write_regions = write_regions

    def __str__(self) -> str:
        return "Processor to perform photometric calibration."

    def get_phot_output_dir(self):
        """
        Return the
        :return:
        """
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def calculate_zeropoint(self, ref_cat_path: Path, img_cat_path: Path) -> list[dict]:
        """
        Function to calculate zero point from two catalogs
        Args:
            ref_cat_path: Path to reference ldac catalog
            img_cat_path: Path to image ldac catalog
        Returns:
        """
        ref_cat = get_table_from_ldac(ref_cat_path)
        img_cat = get_table_from_ldac(img_cat_path)

        if len(ref_cat) == 0:
            err = "No sources found in reference catalog"
            logger.error(err)
            raise PhotometryReferenceError(err)

        ref_coords = SkyCoord(ra=ref_cat["ra"], dec=ref_cat["dec"], unit=(u.deg, u.deg))
        clean_mask = (
            (img_cat["FLAGS"] == 0)
            & (img_cat["FWHM_WORLD"] < self.fwhm_threshold_arcsec / 3600.0)
            & (img_cat["X_IMAGE"] > self.x_lower_limit)
            & (img_cat["X_IMAGE"] < self.x_upper_limit)
            & (img_cat["Y_IMAGE"] > self.y_lower_limit)
            & (img_cat["Y_IMAGE"] < self.y_upper_limit)
        )

        img_coords = SkyCoord(
            ra=img_cat["ALPHAWIN_J2000"],
            dec=img_cat["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )
        clean_img_cat = img_cat[clean_mask]
        logger.debug(f"Found {len(clean_img_cat)} clean sources in image.")
        clean_img_coords = SkyCoord(
            ra=clean_img_cat["ALPHAWIN_J2000"],
            dec=clean_img_cat["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )

        if len(clean_img_coords) == 0:
            err = "No clean sources found in image"
            logger.error(err)
            raise PhotometrySourceError(err)

        idx, d2d, _ = ref_coords.match_to_catalog_sky(clean_img_coords)
        match_mask = d2d < 1.0 * u.arcsec
        matched_ref_cat = ref_cat[match_mask]
        matched_img_cat = clean_img_cat[idx[match_mask]]
        logger.info(
            f"Cross-matched {len(matched_img_cat)} sources from catalog to the image."
        )

        if self.write_regions:
            ref_regions_path = get_output_path(
                base_name="ref.reg",
                dir_root=self.temp_output_sub_dir,
                sub_dir=self.night_sub_dir,
            )
            cleaned_img_regions_path = get_output_path(
                base_name="cleaned_img.reg",
                dir_root=self.temp_output_sub_dir,
                sub_dir=self.night_sub_dir,
            )
            img_regions_path = get_output_path(
                base_name="img.reg",
                dir_root=self.temp_output_sub_dir,
                sub_dir=self.night_sub_dir,
            )

            write_regions_file(
                regions_path=ref_regions_path,
                x_coords=ref_coords.ra.deg,
                y_coords=ref_coords.dec.deg,
                system="wcs",
                region_radius=2.0 / 3600,
            )
            write_regions_file(
                regions_path=cleaned_img_regions_path,
                x_coords=clean_img_coords.ra.deg,
                y_coords=clean_img_coords.dec.deg,
                system="wcs",
                region_radius=2.0 / 3600,
            )
            write_regions_file(
                regions_path=img_regions_path,
                x_coords=img_coords.ra.deg,
                y_coords=img_coords.dec.deg,
                system="wcs",
                region_radius=2.0 / 3600,
            )

        if len(matched_img_cat) < self.num_matches_threshold:
            err = (
                "Not enough cross-matched sources "
                "found to calculate a reliable zeropoint."
            )
            logger.error(err)
            raise PhotometryCrossMatchError(err)

        apertures = self.get_sextractor_apertures()  # aperture diameters
        zeropoints = []

        for i, aperture in enumerate(apertures):
            offsets = np.ma.array(
                matched_ref_cat["magnitude"] - matched_img_cat["MAG_APER"][:, i]
            )
            cl_offset = sigma_clip(offsets)
            num_stars = np.sum(np.invert(cl_offset.mask))

            zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets)

            check = [np.isnan(x) for x in [zp_mean, zp_med, zp_std]]
            if np.sum(check) > 0:
                err = (
                    f"Error with nan when calculating sigma stats: \n "
                    f"mean: {zp_mean}, median: {zp_med}, std: {zp_std}"
                )
                logger.error(err)
                raise PhotometryCalculationError(err)

            zero_dict = {
                "diameter": aperture,
                "zp_mean": zp_mean,
                "zp_median": zp_med,
                "zp_std": zp_std,
                "nstars": num_stars,
                "mag_cat": matched_ref_cat["magnitude"][np.invert(cl_offset.mask)],
                "mag_apers": matched_img_cat["MAG_APER"][:, i][
                    np.invert(cl_offset.mask)
                ],
            }
            zeropoints.append(zero_dict)

        offsets = np.ma.array(
            matched_ref_cat["magnitude"] - matched_img_cat["MAG_AUTO"]
        )
        cl_offset = sigma_clip(offsets, sigma=3)
        num_stars = np.sum(np.invert(cl_offset.mask))
        zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets, sigma=3)
        zero_auto_mag_cat = matched_ref_cat["magnitude"][np.invert(cl_offset.mask)]
        zero_auto_mag_img = matched_img_cat["MAG_AUTO"][np.invert(cl_offset.mask)]
        zeropoints.append(
            {
                "diameter": "AUTO",
                "zp_mean": zp_mean,
                "zp_median": zp_med,
                "zp_std": zp_std,
                "nstars": num_stars,
                "mag_cat": zero_auto_mag_cat,
                "mag_apers": zero_auto_mag_img,
            }
        )

        return zeropoints

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        phot_output_dir = self.get_phot_output_dir()

        try:
            os.makedirs(phot_output_dir)
        except OSError:
            pass

        for image in batch:
            ref_catalog = self.ref_catalog_generator(image)
            ref_cat_path = ref_catalog.write_catalog(image, output_dir=phot_output_dir)
            temp_cat_path = copy_temp_file(
                output_dir=phot_output_dir, file_path=image[SEXTRACTOR_HEADER_KEY]
            )

            temp_files = [temp_cat_path]
            fwhm_med, _, fwhm_std, med_fwhm_pix, _, _ = self.get_fwhm(temp_cat_path)
            image["FWHM_MED"] = fwhm_med
            image["FWHM_STD"] = fwhm_std

            zp_dicts = self.calculate_zeropoint(ref_cat_path, temp_cat_path)

            aperture_diameters = []
            zp_values = []
            for zpvals in zp_dicts:
                image[f"ZP_{zpvals['diameter']}"] = zpvals["zp_mean"]
                image[f"ZP_{zpvals['diameter']}_std"] = zpvals["zp_std"]
                image[f"ZP_{zpvals['diameter']}_nstars"] = zpvals["nstars"]
                try:
                    aperture_diameters.append(float(zpvals["diameter"]))
                    zp_values.append(zpvals["zp_mean"])
                except ValueError:
                    continue

            aperture_diameters.append(med_fwhm_pix)
            zp_values.append(image["ZP_AUTO"])

            if sextractor_checkimg_map["BACKGROUND_RMS"] in image:
                limmags = self.get_maglim(
                    image[sextractor_checkimg_map["BACKGROUND_RMS"]],
                    zp_values,
                    np.array(aperture_diameters) / 2.0,
                )
            else:
                limmags = [-99] * len(aperture_diameters)

            for ind, diam in enumerate(aperture_diameters):
                image[f"MAGLIM_{int(diam)}"] = limmags[ind]
            image["MAGLIM"] = limmags[-1]

            if not self.cache:
                for temp_file in temp_files:
                    temp_file.unlink()
                    logger.debug(f"Deleted temporary file {temp_file}")

        return batch

    @staticmethod
    def get_fwhm(img_cat_path):
        """
        Calculate median FWHM from a ldac path
        Args:
            img_cat_path:

        Returns:
        """
        imcat = get_table_from_ldac(img_cat_path)
        # TODO: de-hardcode
        nemask = (
            (imcat["X_IMAGE"] > 50)
            & (imcat["X_IMAGE"] < 2000)
            & (imcat["Y_IMAGE"] > 50)
            & (imcat["Y_IMAGE"] < 2000)
        )
        imcat = imcat[nemask]
        med_fwhm = np.median(imcat["FWHM_WORLD"])
        mean_fwhm = np.mean(imcat["FWHM_WORLD"])
        std_fwhm = np.std(imcat["FWHM_WORLD"])

        med_fwhm_pix = np.median(imcat["FWHM_IMAGE"])
        mean_fwhm_pix = np.mean(imcat["FWHM_IMAGE"])
        std_fwhm_pix = np.std(imcat["FWHM_IMAGE"])
        return med_fwhm, mean_fwhm, std_fwhm, med_fwhm_pix, mean_fwhm_pix, std_fwhm_pix

    @staticmethod
    def get_maglim(
        bkg_rms_image_path: str | Path,
        zeropoint: float | list[float],
        aperture_radius_pixels: float | list[float],
    ) -> float:
        """
        Function to calculate limiting magnitude
        Args:
            bkg_rms_image_path:
            zeropoint:
            aperture_radius_pixels:
        Returns:
        """
        if isinstance(zeropoint, float):
            zeropoint = [zeropoint]
        if isinstance(aperture_radius_pixels, float):
            aperture_radius_pixels = [aperture_radius_pixels]

        zeropoint = np.array(zeropoint, dtype=float)
        aperture_radius_pixels = np.array(aperture_radius_pixels, dtype=float)
        bkg_rms_image = fits.getdata(bkg_rms_image_path)
        bkg_rms_med = np.nanmedian(bkg_rms_image)
        noise = bkg_rms_med * np.sqrt(np.pi * aperture_radius_pixels)
        maglim = -2.5 * np.log10(5 * noise) + zeropoint
        return maglim

    def get_sextractor_module(self) -> Sextractor:
        mask = [isinstance(x, Sextractor) for x in self.preceding_steps]
        return np.array(self.preceding_steps)[mask][-1]

    def check_prerequisites(
        self,
    ):
        mask = [isinstance(x, Sextractor) for x in self.preceding_steps]
        if np.sum(mask) < 1:
            err = (
                f"{self.__module__} requires {Sextractor} as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise PrerequisiteError(err)

        sextractor_param_path = self.get_sextractor_module().parameters_name

        logger.debug(f"Checking file {sextractor_param_path}")

        with open(sextractor_param_path, "rb") as param_file:
            sextractor_params = [
                x.strip().decode() for x in param_file.readlines() if len(x.strip()) > 0
            ]
            sextractor_params = [
                x.split("(")[0] for x in sextractor_params if x[0] not in ["#"]
            ]

        for param in REQUIRED_PARAMETERS:
            if param not in sextractor_params:
                err = (
                    f"Missing parameter: {self.__module__} requires {param} to run, "
                    f"but this parameter was not found in sextractor config file "
                    f"'{sextractor_param_path}' . "
                    f"Please add the parameter to this list!"
                )
                logger.error(err)
                raise PrerequisiteError(err)

    def get_sextractor_apertures(self) -> list[float]:
        """
        Function to extract sextractor aperture sizes from config file
        Returns:
        """
        sextractor_config_path = self.get_sextractor_module().config

        with open(sextractor_config_path, "rb") as sextractor_config_file:
            aperture_lines = [
                x.decode()
                for x in sextractor_config_file.readlines()
                if np.logical_and(b"PHOT_APERTURES" in x, x.decode()[0] != "#")
            ]

        if len(aperture_lines) > 1:
            err = (
                f"The config file {sextractor_config_path} has "
                f"multiple entries for PHOT_APERTURES."
            )
            logger.error(err)
            raise ProcessorError(err)

        line = aperture_lines[0].replace("PHOT_APERTURES", " ").split("#")[0]

        return [float(x) for x in line.split(",") if x not in [""]]
