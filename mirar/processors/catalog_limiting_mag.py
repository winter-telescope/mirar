"""
Module to calculate the brightest and faintest stars detected in the image
"""

import logging
from pathlib import Path
from typing import Callable

import numpy as np
from astropy.table import Table

from mirar.data import Image, ImageBatch
from mirar.data.utils.coords import write_regions_file
from mirar.paths import SEXTRACTOR_HEADER_KEY, ZP_KEY
from mirar.processors import BaseImageProcessor
from mirar.processors.astromatic.sextractor.sextractor import (
    check_sextractor_prerequisite,
)
from mirar.processors.base_catalog_xmatch_processor import (
    default_image_sextractor_catalog_purifier,
)
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


def default_lim_mag_sextractor_catalog_purifier(
    catalog: Table,
    image: Image,
) -> Table:
    """
    Default function to purify the photometric image catalog
    """
    clean_catalog = default_image_sextractor_catalog_purifier(
        catalog, image, fwhm_threshold_arcsec=20
    )
    for key in ["MAG_AUTO", "SNR_WIN", "FLAGS"]:
        assert key in clean_catalog.keys(), (
            f"Default catalog purifier requires {key},"
            f"which was not found in the catalog. "
            f"Please provide your own "
            f"image_photometric_catalog_purifier"
        )
    clean_catalog_mask = (
        (clean_catalog["MAG_AUTO"] != 99)
        & (clean_catalog["SNR_WIN"] > 0)
        & (clean_catalog["FLAGS"] == 0)
    )
    return clean_catalog[clean_catalog_mask]


class CatalogLimitingMagnitudeCalculator(BaseImageProcessor):
    """
    Processor to calculate the limiting magnitude of an image based on the
    SEXTRACTOR catalog
    """

    base_key = "catlimmagcalc"

    def __init__(
        self,
        image_photometric_catalog_purifier: Callable[
            [Table, Image], Table
        ] = default_lim_mag_sextractor_catalog_purifier,
        sextractor_mag_key_name: str = "MAG_AUTO",
        write_regions: bool = False,
    ):
        super().__init__()
        self.image_photometric_catalog_purifier = image_photometric_catalog_purifier
        self.sextractor_mag_key_name = sextractor_mag_key_name
        self.write_regions = write_regions

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image in batch:
            sextractor_catalog_path = Path(image[SEXTRACTOR_HEADER_KEY])
            image_cat = get_table_from_ldac(sextractor_catalog_path)
            cleaned_image_cat = self.image_photometric_catalog_purifier(
                image_cat, image
            )
            zero_point = image[ZP_KEY]
            detected_mags = cleaned_image_cat[self.sextractor_mag_key_name] + zero_point
            detected_mags_05, detected_mags_50, detected_mags_95 = np.percentile(
                np.array(detected_mags[detected_mags > 0]), [5, 50, 95]
            )
            image["DETMAG95"] = detected_mags_95
            image["DETMAG50"] = detected_mags_50
            image["DETMAG05"] = detected_mags_05

            if self.write_regions:
                ra, dec = (
                    cleaned_image_cat["ALPHAWIN_J2000"],
                    cleaned_image_cat["DELTAWIN_J2000"],
                )
                write_regions_file(
                    x_coords=ra,
                    y_coords=dec,
                    region_radius=4.0 / 3600,
                    regions_path=sextractor_catalog_path.with_suffix(
                        ".catlimmagcalc.reg"
                    ),
                    text=[
                        str(round(x, 2))
                        for x in cleaned_image_cat[self.sextractor_mag_key_name]
                    ],
                )

        return batch

    def check_prerequisites(
        self,
    ):
        check_sextractor_prerequisite(self)
