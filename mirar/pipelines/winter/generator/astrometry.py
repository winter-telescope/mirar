"""
Module for astrometric calibration for the winter pipeline.
"""

import logging
from pathlib import Path

import numpy as np
from astropy.table import Table

from mirar.catalog import CatalogFromFile, Gaia2Mass
from mirar.data import Image
from mirar.paths import REF_CAT_PATH_KEY
from mirar.pipelines.winter.generator.utils import check_winter_local_catalog_overlap

logger = logging.getLogger(__name__)


def winter_astrometry_sextractor_catalog_purifier(catalog: Table, _) -> Table:
    """
    Function to purify the Sextractor catalog for WINTER astrometry
    """
    clean_catalog = catalog[
        (catalog["FLAGS"] == 0) & (catalog["FWHM_IMAGE"] > 0) & (catalog["SNR_WIN"] > 0)
    ]
    return clean_catalog


def winter_astrometric_ref_catalog_generator(
    image: Image,
) -> Gaia2Mass | CatalogFromFile:
    """
    Function to generate a reference catalog for WINTER astrometry

    :param image: Image
    :return: catalogue
    """
    if REF_CAT_PATH_KEY in image.header:
        ref_cat_path = Path(image[REF_CAT_PATH_KEY])
        logger.debug(f"Looking for local reference catalog at {ref_cat_path}")
        if ref_cat_path.exists():
            # Check if there is enough overlap between the image and the
            # local reference catalog
            if check_winter_local_catalog_overlap(ref_cat_path, image):
                logger.debug(f"Loading reference catalog from {ref_cat_path}")
                return CatalogFromFile(catalog_path=ref_cat_path)

            logger.debug(
                "More than 50% of the local reference catalog is "
                "outside the image. Requerying."
            )

    search_radius_arcmin = (
        np.max([image["NAXIS1"], image["NAXIS2"]])
        * np.max([np.abs(image["CD1_1"]), np.abs(image["CD1_2"])])
        * 60
    ) / 2.0
    return Gaia2Mass(
        min_mag=7,
        max_mag=20,
        search_radius_arcmin=search_radius_arcmin,
        cache_catalog_locally=True,
    )
