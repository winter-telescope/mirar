# from mirar.data import Image
import logging

import numpy as np
from astropy.table import Table

from mirar.catalog import PS1, Gaia2Mass
from mirar.data import Image
from mirar.paths import FILTER_KEY
from mirar.pipelines.spring.config import sextractor_astrometry_config
from mirar.processors.base_catalog_xmatch_processor import (
    default_image_sextractor_catalog_purifier,
)

logger = logging.getLogger(__name__)


def spring_anet_sextractor_config_path_generator(_image: Image) -> str:
    """
    Generate the path to the ANET SExtractor configuration file for SPRING images
    Parameters
    ----------
    image:Image

    Returns
    -------
    sextractor_config_path:str

    """
    return sextractor_astrometry_config["config_path"]


def spring_photometric_catalog_generator(image: Image) -> Gaia2Mass | PS1:
    """
    Function to match SPRING image to GAIA/2MASS/PS1 for photometry
    Parameters
    ----------
    image:Image

    Returns
    ----------
    catalogue
    """
    filter_name = image[FILTER_KEY]
    search_radius_arcmin = (
        np.max([image["NAXIS1"], image["NAXIS2"]])
        * np.max([np.abs(image["CD1_1"]), np.abs(image["CD1_2"])])
        * 60.0
    ) / 2.0

    if filter_name in ["J", "H"]:
        return Gaia2Mass(
            min_mag=0,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name,
            snr_threshold=20,
            cache_catalog_locally=False,
        )

    if filter_name in ["Y"]:
        return PS1(
            min_mag=0,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name.lower(),
            cache_catalog_locally=False,
        )

    err = f"Filter name {filter_name} not recognized"
    logger.error(err)
    raise ValueError(err)


def spring_ref_photometric_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> tuple[Table, Table]:
    """
    Default function to purify the photometric image catalog
    """
    return default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog,
        ref_catalog=ref_catalog,
        image=image,
        edge_width_pixels=100,
        fwhm_threshold_arcsec=6,
    )


def spring_photcal_color_columns_generator(image: Image) -> tuple[list, list, tuple]:
    filter_name = image[FILTER_KEY]
    if filter_name in ["Y"]:
        return ["ymag", "zmag"], ["e_ymag", "e_zmag"], (0, 25)
    if filter_name == "J":
        return ["j_m", "h_m"], ["j_msigcom", "h_msigcom"], (0, 25)
    if filter_name == "H":
        return ["h_m", "k_m"], ["h_msigcom", "k_msigcom"], (0, 25)
    err = f"Filter name {filter_name} not recognized"
    logger.error(err)
    raise ValueError(err)
