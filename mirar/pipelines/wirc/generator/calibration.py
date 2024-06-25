"""
Module containing WIRC-specific generator functions to
yield e.g catalog for astrometric calibrations
"""

import logging

from mirar.catalog import Gaia2Mass
from mirar.data import Image

logger = logging.getLogger(__name__)


def wirc_astrometric_catalog_generator(_) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for astrometry

    :return: catalogue
    """
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=10)


def wirc_photometric_catalog_generator(image: Image) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    filter_name = image["FILTER"]
    if filter_name == "ks":
        return Gaia2Mass(
            min_mag=10, max_mag=20, search_radius_arcmin=10, filter_name="k"
        )
    return Gaia2Mass(
        min_mag=10, max_mag=20, search_radius_arcmin=10, filter_name=filter_name
    )
