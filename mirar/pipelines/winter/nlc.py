"""
Module for applying WINTER non-linear correction to images
"""

from winternlc.non_linear_correction import apply_nlc_single

from mirar.data import ImageBatch


def apply_winter_nlc(images: ImageBatch) -> ImageBatch:
    """
    Apply WINTER non-linear correction to images.
    Uses header information to apply the correct correction.

    :param images: ImageBatch to apply non-linear correction to
    :return: Corrected ImageBatch
    """

    for image in images:
        data = image.get_data()
        header = image.get_header()
        corrected_data = apply_nlc_single(data, header)
        image.set_data(corrected_data)

    return images
