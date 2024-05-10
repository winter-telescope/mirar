"""
Module for photometric for the Winter pipeline.
"""

import logging
from pathlib import Path

import numpy as np
from astropy.table import Table

from mirar.catalog import PS1, CatalogFromFile, Gaia2Mass
from mirar.data import Image
from mirar.paths import FILTER_KEY, REF_CAT_PATH_KEY
from mirar.pipelines.winter.generator.utils import check_winter_local_catalog_overlap
from mirar.processors.base_catalog_xmatch_processor import (
    default_image_sextractor_catalog_purifier,
)
from mirar.processors.photcal.photcalibrator import PhotCalibrator

logger = logging.getLogger(__name__)


def winter_astrostat_catalog_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> (Table, Table):
    """
    Default function to purify the photometric image catalog
    """

    return default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog,
        ref_catalog=ref_catalog,
        image=image,
        edge_width_pixels=0,
        fwhm_threshold_arcsec=20.0,
    )


def winter_photometric_catalog_generator(
    image: Image,
) -> Gaia2Mass | PS1 | CatalogFromFile:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    if REF_CAT_PATH_KEY in image.header:
        ref_cat_path = Path(image[REF_CAT_PATH_KEY])
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

    filter_name = image["FILTER"]
    search_radius_arcmin = (
        np.max([image["NAXIS1"], image["NAXIS2"]])
        * np.max([np.abs(image["CD1_1"]), np.abs(image["CD1_2"])])
        * 60
    ) / 2.0

    if filter_name in ["J", "H"]:
        return Gaia2Mass(
            min_mag=0,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name,
            snr_threshold=20,
            cache_catalog_locally=True,
        )

    if filter_name in ["Y"]:
        return PS1(
            min_mag=0,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name.lower(),
            cache_catalog_locally=True,
        )

    err = f"Filter {filter_name} not recognised"
    logger.error(err)
    raise ValueError(err)


def winter_photometric_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> (Table, Table):
    """
    Default function to purify the photometric image catalog
    """

    sci_catalog, ref_catalog = default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog,
        ref_catalog=ref_catalog,
        image=image,
        edge_width_pixels=100,
        fwhm_threshold_arcsec=5.0,
    )

    sci_catalog = sci_catalog[
        (sci_catalog["FLAGS_MODEL"] == 0)
        & (sci_catalog["FLUX_MAX"] < 30000)
        & (sci_catalog["FLUX_MAX"] > 0)
    ]

    ref_catalog = ref_catalog[ref_catalog["magnitude"] > 10]
    return sci_catalog, ref_catalog


def winter_photcal_color_columns_generator(image):
    """
    Returns the color columns for WINTER photometric calibration

    :param image: Image
    :return: color columns
    """
    filter_name = image[FILTER_KEY]
    if filter_name == "J":
        return ["j_m", "h_m"], ["j_msigcom", "h_msigcom"], (0, 25)
    if filter_name == "H":
        return ["h_m", "k_m"], ["h_msigcom", "k_msigcom"], (0, 25)
    if filter_name in ["Y"]:
        return ["ymag", "zmag"], ["e_ymag", "e_zmag"], (0, 25)
    err = f"Filter {filter_name} not recognised"
    logger.error(err)
    raise ValueError(err)


def winter_ref_photometric_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> (Table, Table):
    """
    Default function to purify the photometric image catalog
    """

    return default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog,
        ref_catalog=ref_catalog,
        image=image,
        edge_width_pixels=100,
        fwhm_threshold_arcsec=4.0,
    )


def winter_reference_phot_calibrator(_: Image, **kwargs) -> PhotCalibrator:
    """
    Generates a resampler for reference images

    :param _: image
    :param kwargs: kwargs
    :return: Swarp processor
    """

    return PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        write_regions=True,
        catalogs_purifier=winter_ref_photometric_catalogs_purifier,
        **kwargs,
    )
