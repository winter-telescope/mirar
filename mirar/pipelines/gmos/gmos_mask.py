"""
Functions to create a mask for GMOS data.
"""

import logging

import numpy as np

from mirar.data import Image

logger = logging.getLogger(__name__)


def generate_gmos_mask(image: Image) -> np.array:
    """
    Function to create a mask for GMOS data

    :param image: Image object
    :return: Image object with mask
    """

    data = image.get_data()

    median = np.nanmedian(data)
    std = np.nanstd(data)

    print("Median: ", median)
    print("STD: ", std)

    mask = data < 0.5 * median

    print("Mask: ", mask)

    logger.info(f"Masking {np.sum(mask)}/{len(mask)} pixels")

    return mask
