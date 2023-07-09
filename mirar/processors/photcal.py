"""
Module for running photometric calibration
"""
import logging
import os
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.table import Table

from mirar.catalog.base_catalog import BaseCatalog
from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import get_output_dir
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
    bkg_rms_image = fits.getdata(bkg_rms_image_path)
    bkg_rms_med = np.nanmedian(bkg_rms_image)
    noise = bkg_rms_med * np.sqrt(np.pi * aperture_radius_pixels)
    maglim = -2.5 * np.log10(5 * noise) + zeropoint
    return maglim


class PhotCalibrator(BaseProcessorWithCrossMatch):
    """
    Photometric calibrator processor
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

    def __str__(self) -> str:
        return "Processor to perform photometric calibration."

    def get_phot_output_dir(self):
        """
        Return the
        :return:
        """
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def calculate_zeropoint(
        self,
        ref_cat: Table,
        clean_img_cat: Table,
    ) -> list[dict]:
        """
        Function to calculate zero point from two catalogs
        Args:
            ref_cat: Reference catalog table
            clean_img_cat: Catalog of sources from image to xmatch with ref_cat
        Returns:
        """

        matched_img_cat, matched_ref_cat, _ = self.xmatch_catalogs(
            ref_cat=ref_cat,
            image_cat=clean_img_cat,
            crossmatch_radius_arcsec=self.crossmatch_radius_arcsec,
        )
        logger.debug(
            f"Cross-matched {len(matched_img_cat)} sources from catalog to the image."
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
            ref_cat, _, cleaned_img_cat = self.setup_catalogs(image)

            fwhm_med, _, fwhm_std, med_fwhm_pix, _, _ = get_fwhm(cleaned_img_cat)
            image["FWHM_MED"] = fwhm_med
            image["FWHM_STD"] = fwhm_std
            image["FWHM_PIX"] = med_fwhm_pix

            if len(ref_cat) == 0:
                err = "No sources found in reference catalog"
                logger.error(err)
                raise PhotometryReferenceError(err)

            logger.debug(f"Found {len(cleaned_img_cat)} clean sources in image.")

            if len(cleaned_img_cat) == 0:
                err = "No clean sources found in image"
                logger.error(err)
                raise PhotometrySourceError(err)

            zp_dicts = self.calculate_zeropoint(
                ref_cat=ref_cat, clean_img_cat=cleaned_img_cat
            )

            aperture_diameters = []
            zp_values = []

            with warnings.catch_warnings(record=True):
                warnings.simplefilter("ignore", category=VerifyWarning)

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

                for ind, diam in enumerate(aperture_diameters):
                    image[f"MAGLIM_{int(diam)}"] = limmags[ind]
                image["MAGLIM"] = limmags[-1]

        return batch
