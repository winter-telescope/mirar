"""
Some utility functions for the WINTER pipeline
"""

import logging
from pathlib import Path

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.wcs import NoConvergence

from mirar.data import Image
from mirar.data.utils.coords import check_coords_within_image
from mirar.paths import TIME_KEY
from mirar.pipelines.winter.models import DEFAULT_FIELD
from mirar.processors.flat import get_convolution
from mirar.utils import get_table_from_ldac

logger = logging.getLogger(__name__)


def check_winter_local_catalog_overlap(ref_cat_path: Path, image: Image) -> bool:
    """
    Returns a function that returns the local reference catalog
    """
    # Check if there is enough overlap between the image and the
    # local reference catalog
    local_ref_cat = get_table_from_ldac(ref_cat_path)

    if len(local_ref_cat) == 0:
        logger.debug(f"Reference catalog {ref_cat_path} is empty.")
        return False

    try:
        srcs_in_image = check_coords_within_image(
            ra=local_ref_cat["ra"], dec=local_ref_cat["dec"], header=image.get_header()
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
    return cat_overlaps


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
            f"_{image['TARGNAME']}_{image[TIME_KEY]}"
            f"_{image['FILTER']}.ldac.cat"
        )
    return ref_cat_path


def get_outlier_pixel_mask(img: np.ndarray, thresh: float = 3.0) -> np.ndarray:
    """
    Get oulier pixels that are above or below a threshold
    :param img: np.ndarray
    :param thresh: float
    """
    _, median, std = sigma_clipped_stats(img, sigma=3.0)
    return (img < median - thresh * std) | (img > median + thresh * std)


def construct_smooth_gradient_for_image(data: np.ndarray) -> np.ndarray:
    """
    Construct a smooth gradient for the image
    :param data: np.ndarray
    :return: np.ndarray
    """
    smooth_img = get_convolution(data, 100)
    return smooth_img


def smooth_and_normalize_image(data: np.ndarray) -> np.ndarray:
    """
    Smooth and normalize the image
    :param data: np.ndarray
    :return: np.ndarray
    """
    smooth_img = get_convolution(data, 100)
    return data / smooth_img


def get_smoothened_outlier_pixel_mask_from_list(
    img_data_list: list[np.ndarray],
) -> np.ndarray:
    """
    Take a list of images and return a mask of outlier pixels after removing a smooth
    gradient from them
    :param img_data_list: list[np.ndarray]
    :return: np.ndarray
    """
    masks = []
    for data in img_data_list:
        smooth_img = smooth_and_normalize_image(data)
        mask = get_outlier_pixel_mask(smooth_img, thresh=3.0)
        masks.append(mask)
    return np.logical_and.reduce(masks)
