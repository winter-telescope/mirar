"""
Module for applying WINTER non-linear correction to images
"""

import logging

from winternlc import apply_nlc_single

from mirar.data import ImageBatch

logger = logging.getLogger(__name__)


def apply_winter_nlc(images: ImageBatch) -> ImageBatch:
    """
    Apply WINTER non-linear correction to images.
    Uses header information to apply the correct correction.

    :param images: ImageBatch to apply non-linear correction to
    :return: Corrected ImageBatch
    """
    logger.debug("Applying WINTER non-linear correction")
    for image in images:
        data = image.get_data()
        header = image.get_header()
        if int(header["BOARD_ID"]) in [2, 3, 4, 5, 6]:
            corrected_data = apply_nlc_single(data, header)
            image.set_data(corrected_data)

    return images
