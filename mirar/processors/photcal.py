"""
Module for running photometric calibration
"""
import logging
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.table import Table
from scipy.optimize import curve_fit

from mirar.catalog.base_catalog import BaseCatalog
from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import MAGLIM_KEY, ZP_KEY, ZP_NSTARS_KEY, ZP_STD_KEY, get_output_dir
from mirar.processors.astromatic.sextractor.sextractor import sextractor_checkimg_map
from mirar.processors.astrometry.validate import get_fwhm
from mirar.processors.base_catalog_xmatch_processor import (
    BaseProcessorWithCrossMatch,
    default_image_sextractor_catalog_purifier,
)

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
    logger.debug(aperture_radius_pixels)
    bkg_rms_image = fits.getdata(bkg_rms_image_path)
    bkg_rms_med = np.nanmedian(bkg_rms_image)
    noise = bkg_rms_med * np.sqrt(np.pi * aperture_radius_pixels**2)
    maglim = -2.5 * np.log10(5 * noise) + zeropoint
    logger.debug(f"Aperture radii: {aperture_radius_pixels}")
    logger.debug(f"Calculated maglim: {maglim}")
    return maglim


class BaseZeroPointCalculator:
    """
    Base class for actually calculating zero point calculators
    """

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        """
        Function to calculate zero point from two catalogs and add the zeropoint
        information to the image header
        Args:
            matched_ref_cat: Reference catalog table
            matched_img_cat: Catalog of sources
            image: Image object
            colnames: List of column names from the image catalog to use for
            calculating zero point. The reference catalog is assumed to have the
            "magnitude" column, as it comes from mirar.catalog.base_catalog.BaseCatalog
        Returns:
            Image object with zero point information added to the header.
        """
        raise NotImplementedError


class OutlierRejectionZPCalculator(BaseZeroPointCalculator):
    """
    Class to calculate zero point using outlier rejection
    Attributes:
        num_stars_threshold: int to use as minimum number of stars
        outlier_rejection_threshold: float or list of floats to use as number of sigmas
    """

    def __init__(
        self,
        num_stars_threshold: int = 5,
        outlier_rejection_threshold: float | list[float] = 3.0,
    ):
        self.num_stars_threshold = num_stars_threshold
        self.outlier_rejection_threshold = outlier_rejection_threshold
        if isinstance(outlier_rejection_threshold, float):
            self.outlier_rejection_threshold = [outlier_rejection_threshold]
        self.outlier_rejection_threshold = np.sort(self.outlier_rejection_threshold)

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        zeropoints = []

        for colname in colnames:
            offsets = np.ma.array(
                matched_ref_cat["magnitude"] - matched_img_cat[colname]
            )
            for outlier_thresh in self.outlier_rejection_threshold:
                cl_offset = sigma_clip(offsets, sigma=outlier_thresh)
                num_stars = np.sum(np.invert(cl_offset.mask))

                zp_mean, zp_med, zp_std = sigma_clipped_stats(
                    offsets, sigma=outlier_thresh
                )

                if num_stars > self.num_stars_threshold:
                    break

            check = [np.isnan(x) for x in [zp_mean, zp_med, zp_std]]
            if np.sum(check) > 0:
                err = (
                    f"Error with nan when calculating sigma stats: \n "
                    f"mean: {zp_mean}, median: {zp_med}, std: {zp_std}"
                )
                logger.error(err)
                raise PhotometryCalculationError(err)

            aperture = colname.split("_")[-1]
            zero_dict = {
                "diameter": aperture,
                "zp_mean": zp_mean,
                "zp_median": zp_med,
                "zp_std": zp_std,
                "nstars": num_stars,
                "mag_cat": matched_ref_cat["magnitude"][np.invert(cl_offset.mask)],
                "mag_apers": matched_img_cat[colname][np.invert(cl_offset.mask)],
            }
            zeropoints.append(zero_dict)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore", category=VerifyWarning)

            for zpvals in zeropoints:
                image[f"ZP_{zpvals['diameter']}"] = zpvals["zp_mean"]
                image[f"ZP_{zpvals['diameter']}_std"] = zpvals["zp_std"]
                image[f"ZP_{zpvals['diameter']}_nstars"] = zpvals["nstars"]

        return image


def line_model(data, slope, intercept):
    """linear model to hand to scipy curve_fit"""
    return slope * data + intercept


class ZPWithColorTermCalculator(BaseZeroPointCalculator):
    """
    Class to calculate zero point by including a color term. This models the data as

    ref_mag - img_mag = ZP + C * (ref_color)

    Attributes:
        color_colnames_generator: function that takes an image as input and returns
        a list containing two strings that are the column names of the reference catalog
        to use for the color term. The first string is the bluer band, the second is the
         redder band.
    """

    def __init__(
        self,
        color_colnames_generator: Callable[Image, [list[str, str]]],
    ):
        self.color_colnames_generator = color_colnames_generator

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        color_colnames = self.color_colnames_generator(image)
        colors = matched_ref_cat[color_colnames[0]] - matched_ref_cat[color_colnames[1]]

        for colname in colnames:
            y = matched_img_cat[colname] - matched_ref_cat["magnitude"]
            x = colors
            y_err = np.sqrt(
                matched_img_cat[colname.replace("MAG", "MAGERR")] ** 2
                + matched_ref_cat["magnitude_err"] ** 2
            )
            popt, pcov = curve_fit(f=line_model, xdata=x, ydata=y, sigma=y_err)
            color, zero_point = popt
            color_err, zp_err = np.sqrt(np.diag(pcov))

            aperture = colname.split("_")[-1]
            image[f"ZP_{aperture}"] = zero_point
            image[f"ZP_{aperture}_std"] = zp_err
            image[f"ZP_{aperture}_nstars"] = len(matched_ref_cat)
            image[f"C_{aperture}"] = color
            image[f"C_{aperture}_std"] = color_err

        return image


class PhotCalibrator(BaseProcessorWithCrossMatch):
    """
    Photometric calibrator processor

    Attributes:
        num_matches_threshold: minimum number of matches required for
        photometric calibration
        for outlier rejection. If a ist is provided, the list is sorted and stepped
        through in order with increasing thresholds until the specified
        number of matches is reached.
    """

    base_key = "photcalibrator"

    def __init__(
        self,
        ref_catalog_generator: Callable[[Image], BaseCatalog],
        temp_output_sub_dir: str = "phot",
        image_photometric_catalog_purifier: Callable[
            [Table, Image], Table
        ] = default_image_sextractor_catalog_purifier,
        num_matches_threshold: int = 5,
        crossmatch_radius_arcsec: float = 1.0,
        write_regions: bool = False,
        cache: bool = False,
        zp_calculator: BaseZeroPointCalculator = OutlierRejectionZPCalculator(),
    ):
        super().__init__(
            ref_catalog_generator=ref_catalog_generator,
            temp_output_sub_dir=temp_output_sub_dir,
            crossmatch_radius_arcsec=crossmatch_radius_arcsec,
            sextractor_catalog_purifier=image_photometric_catalog_purifier,
            write_regions=write_regions,
            cache=cache,
            required_parameters=REQUIRED_PARAMETERS,
        )
        self.num_matches_threshold = num_matches_threshold
        self.zp_calculator = zp_calculator

    def __str__(self) -> str:
        return "Processor to perform photometric calibration."

    def get_phot_output_dir(self):
        """
        Return the
        :return:
        """
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        phot_output_dir = self.get_phot_output_dir()
        phot_output_dir.mkdir(parents=True, exist_ok=True)

        apertures = self.get_sextractor_apertures()  # aperture diameters

        for image in batch:
            ref_cat, _, cleaned_img_cat = self.setup_catalogs(image)

            fwhm_med, _, fwhm_std, med_fwhm_pix, _, _ = get_fwhm(cleaned_img_cat)

            header_map = {
                "FWHM_MED": fwhm_med,
                "FWHM_STD": fwhm_std,
                "FWHM_PIX": med_fwhm_pix,
            }
            for key, value in header_map.items():
                if np.isnan(value):
                    value = -999.0
                image.header[key] = value

            if len(ref_cat) < self.num_matches_threshold:
                err = (
                    f"Not enough sources ({len(ref_cat)} found in reference catalog "
                    f"to calculate a reliable zeropoint. "
                    f"Require at least {self.num_matches_threshold} matches."
                )
                logger.error(err)
                raise PhotometryReferenceError(err)

            logger.debug(f"Found {len(cleaned_img_cat)} clean sources in image.")

            if len(cleaned_img_cat) < self.num_matches_threshold:
                err = (
                    f"Not enough sources ({len(cleaned_img_cat)} "
                    f"found in source catalog "
                    f"to calculate a reliable zeropoint. "
                    f"Require at least {self.num_matches_threshold} matches."
                )
                logger.error(err)
                raise PhotometrySourceError(err)

            matched_img_cat, matched_ref_cat, _ = self.xmatch_catalogs(
                ref_cat=ref_cat,
                image_cat=cleaned_img_cat,
                crossmatch_radius_arcsec=self.crossmatch_radius_arcsec,
            )
            logger.debug(
                f"Cross-matched {len(matched_img_cat)} sources from catalog to "
                "the image."
            )

            if len(matched_img_cat) < self.num_matches_threshold:
                err = (
                    "Not enough cross-matched sources "
                    "found to calculate a reliable zeropoint. "
                    f"Only found {len(matched_img_cat)} crossmatches, "
                    f"while {self.num_matches_threshold} are required. "
                    f"Used {len(ref_cat)} reference sources and "
                    f"{len(cleaned_img_cat)} image sources."
                )
                logger.error(err)
                raise PhotometryCrossMatchError(err)

            # Add columns to image catalog for each aperture
            colnames = []
            for ind, aperture in enumerate(apertures):
                matched_img_cat[f"MAGAPER_{aperture}"] = matched_img_cat["MAG_APER"][
                    :, ind
                ]
                matched_img_cat[f"MAGAPERERR_{aperture}"] = matched_img_cat[
                    "MAGERR_APER"
                ][:, ind]

                colnames.append(f"MAGAPER_{aperture}")

            if "MAG_AUTO" in matched_img_cat.colnames:
                colnames.append("MAG_AUTO")

            if "MAG_PSF" in matched_img_cat.colnames:
                colnames.append("MAG_PSF")

            image = self.zp_calculator.calculate_zeropoint(
                image=image,
                matched_ref_cat=matched_ref_cat,
                matched_img_cat=matched_img_cat,
                colnames=colnames,
            )

            aperture_diameters = []
            zp_values = []

            with warnings.catch_warnings(record=True):
                warnings.simplefilter("ignore", category=VerifyWarning)

                for col in colnames:
                    aper = col.split("_")[-1]
                    # Check if the right zeropoint keys are in the image header
                    for key in [f"ZP_{aper}", f"ZP_{aper}_std", f"ZP_{aper}_nstars"]:
                        assert (
                            key in image.header.keys()
                        ), f"Zeropoint key {key} not found in image header."
                    zp_values.append(image[f"ZP_{aper}"])
                    if col in ["MAG_AUTO", "MAG_PSF"]:
                        aperture_diameters.append(med_fwhm_pix * 2)
                    else:
                        aperture_diameters.append(float(aper))
                if sextractor_checkimg_map["BACKGROUND_RMS"] in image.header.keys():
                    logger.debug(
                        "Calculating limiting magnitudes from background RMS file"
                    )
                    limmags = get_maglim(
                        image[sextractor_checkimg_map["BACKGROUND_RMS"]],
                        zp_values,
                        np.array(aperture_diameters) / 2.0,
                    )
                else:
                    limmags = [-99] * len(aperture_diameters)

                for ind, diam in enumerate(aperture_diameters[:-1]):
                    image[f"MAGLIM_{np.rint(diam)}"] = limmags[ind]

                image[MAGLIM_KEY] = limmags[-1]

                zp_key = "AUTO"
                if "MAG_PSF" in colnames:
                    zp_key = "PSF"
                image[ZP_KEY] = image[f"ZP_{zp_key}"]
                image[ZP_STD_KEY] = image[f"ZP_{zp_key}_STD"]
                image[ZP_NSTARS_KEY] = image[f"ZP_{zp_key}_NSTARS"]
                image["MAGSYS"] = "AB"

        return batch
