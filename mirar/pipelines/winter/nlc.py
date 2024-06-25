from winternlc.non_linear_correction import nlc_single

from mirar.data import ImageBatch


def apply_winter_nlc(images: ImageBatch) -> ImageBatch:
    """
    Apply WINTER non-linear correction to images

    :param images: ImageBatch to apply non-linear correction to
    :return: Corrected ImageBatch
    """

    for image in images:
        data = image.get_data()
        board_id = image["BOARD_ID"]
        corrected_image = nlc_single(data, board_id)
        image.set_data(corrected_image)

    return images
