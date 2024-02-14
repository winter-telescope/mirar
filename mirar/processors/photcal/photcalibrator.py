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
from astropy.table import Table

from mirar.catalog.base_catalog import BaseCatalog
from mirar.data import Image, ImageBatch
from mirar.paths import MAGLIM_KEY, ZP_KEY, ZP_NSTARS_KEY, ZP_STD_KEY, get_output_dir
from mirar.processors.astromatic.sextractor.sextractor import sextractor_checkimg_map
from mirar.processors.astrometry.validate import get_fwhm
from mirar.processors.base_catalog_xmatch_processor import (
    BaseProcessorWithCrossMatch,
    default_image_sextractor_catalog_purifier,
)
from mirar.processors.photcal.photcal_errors import (
    PhotometryCrossMatchError,
    PhotometryReferenceError,
    PhotometrySourceError,
)
from mirar.processors.photcal.zp_calculator.base_zp_calculator import (
    BaseZeroPointCalculator,
)
from mirar.processors.photcal.zp_calculator.outlier_rejection_zp_calculator import (
    OutlierRejectionZPCalculator,
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
                matched_img_cat[f"MAGERRAPER_{aperture}"] = matched_img_cat[
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
