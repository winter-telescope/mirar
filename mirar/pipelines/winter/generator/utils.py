"""
Some utility functions for the WINTER pipeline
"""

import logging
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import NoConvergence

from mirar.data import Image
from mirar.data.utils.coords import (
    check_coords_within_image,
    get_corners_ra_dec_from_header,
)
from mirar.paths import TARGET_KEY
from mirar.pipelines.winter.constants import WINTER_DITHER_RADIUS_ARCMIN
from mirar.pipelines.winter.models import DEFAULT_FIELD
from mirar.utils import get_table_from_ldac

logger = logging.getLogger(__name__)


def check_winter_local_catalog_overlap(ref_cat_path: Path, image: Image) -> bool:
    """
    Checks whether a locally-cached reference catalog has adequate coverage
    of `image` to be safely reused, instead of re-querying.

    A cached catalog is keyed by field/subdetector/filter and gets reused
    for every dither of a visit, not just the one image it was originally
    queried for. Checking only that most of the cached stars fall somewhere
    within this image (the bulk-density check below) isn't enough on its
    own: a catalog queried with just enough radius for one dither's own
    footprint can still pass that check while leaving an edge of a
    different, offset dither with no coverage at all - exactly the failure
    mode this also checks for directly, by requiring a nearby reference
    star near each corner of the image.
    """
    local_ref_cat = get_table_from_ldac(ref_cat_path)

    if len(local_ref_cat) == 0:
        logger.debug(f"Reference catalog {ref_cat_path} is empty.")
        return False

    header = image.get_header()

    try:
        srcs_in_image = check_coords_within_image(
            ra=local_ref_cat["ra"], dec=local_ref_cat["dec"], header=header
        )
    except NoConvergence:
        logger.debug(
            f"Reference catalog {ref_cat_path} does not overlap with image."
            "It is so off, that WCS failed to converge on a solution for "
            "the pixels in the image."
        )
        return False

    num_srcs_in_image = np.sum(srcs_in_image)

    cat_overlaps = num_srcs_in_image > len(local_ref_cat) * 0.5
    if not cat_overlaps:
        logger.debug(
            "More than 50% of the local reference catalog is outside the image."
        )
        return False

    cat_coords = SkyCoord(
        ra=local_ref_cat["ra"] * u.deg, dec=local_ref_cat["dec"] * u.deg
    )
    for corner_ra, corner_dec in get_corners_ra_dec_from_header(header):
        corner_coord = SkyCoord(ra=corner_ra * u.deg, dec=corner_dec * u.deg)
        if (
            corner_coord.separation(cat_coords).min()
            > WINTER_DITHER_RADIUS_ARCMIN * u.arcmin
        ):
            logger.debug(
                f"Local reference catalog {ref_cat_path} has no reference "
                f"star within {WINTER_DITHER_RADIUS_ARCMIN} arcmin of "
                f"corner (ra={corner_ra:.4f}, dec={corner_dec:.4f}) of the "
                "image - insufficient coverage, requerying."
            )
            return False

    return True


def winter_ref_catalog_namer(image: Image, output_dir: Path) -> Path:
    """
    Function to name the reference catalog to use for WINTER astrometry

    :param image: Image
    :param output_dir: Path
    :return: Output path
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    if image["FIELDID"] != DEFAULT_FIELD:
        ref_cat_path = (
            output_dir / f"field{image['FIELDID']}_{image['SUBDETID']}"
            f"_{image['FILTER']}.ldac.cat"
        )
    else:
        ref_cat_path = (
            output_dir / f"field{image['FIELDID']}_{image['SUBDETID']}_"
            f"_{image['TARGNAME']}_{image[TARGET_KEY]}"
            f"_{image['FILTER']}.ldac.cat"
        )
    return ref_cat_path
