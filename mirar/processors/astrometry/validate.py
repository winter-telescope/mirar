"""
Module for validating astrometric solutions
"""
import logging
from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.table import Table

from mirar.catalog.base_catalog import BaseCatalog
from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import get_output_dir
from mirar.processors.base_catalog_xmatch_processor import (
    BaseProcessorWithCrossMatch,
    default_image_sextractor_catalog_purifier,
)


class AstrometryValidateCrossmatchError(ProcessorError):
    """
    Class for errors in AstrometryValidate
    """


logger = logging.getLogger(__name__)

# All the Sextractor parameters required for this script to run
REQUIRED_PARAMETERS = [
    "X_IMAGE",
    "Y_IMAGE",
    "FWHM_WORLD",
    "FWHM_IMAGE",
    "FLAGS",
    "ALPHAWIN_J2000",
    "DELTAWIN_J2000",
]


def get_fwhm(cleaned_img_cat: Table):
    """
    Calculate median FWHM from a ldac path
    Args:


    Returns:
    """
    mean_fwhm, med_fwhm, std_fwhm = sigma_clipped_stats(cleaned_img_cat["FWHM_WORLD"])

    mean_fwhm_pix, med_fwhm_pix, std_fwhm_pix = sigma_clipped_stats(
        cleaned_img_cat["FWHM_IMAGE"]
    )
    return med_fwhm, mean_fwhm, std_fwhm, med_fwhm_pix, mean_fwhm_pix, std_fwhm_pix


class AstrometryStatsWriter(BaseProcessorWithCrossMatch):
    """
    Processor to calculate astrometry statistics
    """

    def __init__(
        self,
        ref_catalog_generator: Callable[[Image], BaseCatalog],
        temp_output_sub_dir: str = "astrstat",
        image_catalog_purifier: Callable[
            [Table, Image], Table
        ] = default_image_sextractor_catalog_purifier,
        crossmatch_radius_arcsec: float = 3.0,
        write_regions: bool = False,
        cache: bool = False,
    ):
        super().__init__(
            ref_catalog_generator=ref_catalog_generator,
            temp_output_sub_dir=temp_output_sub_dir,
            crossmatch_radius_arcsec=crossmatch_radius_arcsec,
            sextractor_catalog_purifier=image_catalog_purifier,
            write_regions=write_regions,
            cache=cache,
            required_parameters=REQUIRED_PARAMETERS,
        )

    def __str__(self):
        return (
            "Processor to calculate astrometry statistics and "
            "update the image headers"
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        output_dir = get_output_dir(
            dir_root=self.temp_output_sub_dir, sub_dir=self.night_sub_dir
        )
        if not Path(output_dir).exists():
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        for image in batch:
            ref_cat, _, cleaned_img_cat = self.setup_catalogs(image)
            if len(cleaned_img_cat) == 0:
                warn = "No sources found after cleaning catalog"
                logger.warning(warn)
                fwhm_med, fwhm_std, med_fwhm_pix = -999, -999, -999
            else:
                fwhm_med, _, fwhm_std, med_fwhm_pix, _, _ = get_fwhm(cleaned_img_cat)

                if np.isnan(fwhm_med):
                    fwhm_med = -999
                if np.isnan(fwhm_std):
                    fwhm_std = -999
                if np.isnan(med_fwhm_pix):
                    med_fwhm_pix = -999

            logger.debug(f"FWHM_MED: {fwhm_med}, {len(cleaned_img_cat)}")
            image["FWHM_MED"] = fwhm_med
            image["FWHM_STD"] = fwhm_std
            image["FWHM_PIX"] = med_fwhm_pix

            image.header["ASTUNC"] = -999.0
            image.header["ASTFIELD"] = -999.0

            if (len(ref_cat) > 0) & (len(cleaned_img_cat) > 0):
                _, _, d2d = self.xmatch_catalogs(
                    ref_cat=ref_cat,
                    image_cat=cleaned_img_cat,
                    crossmatch_radius_arcsec=self.crossmatch_radius_arcsec,
                )

                if len(d2d) > 0:
                    image.header["ASTUNC"] = np.nanmedian(d2d.value)
                    image.header["ASTFIELD"] = np.arctan(
                        image.header["CD1_2"] / image.header["CD1_1"]
                    ) * (180 / np.pi)

        return batch
