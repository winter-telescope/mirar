"""
Module to annotate target coordinates on images.
"""

from mirar.data import ImageBatch
from mirar.paths import TARGET_KEY, TIME_KEY


def annotate_target_coordinates(image_batch: ImageBatch) -> ImageBatch:
    """
    Function to annotate target coordinates on images.
    For WASP, this should be the value of RA/DEC in the header of the first image.

    :param image_batch: ImageBatch object
    :return: ImageBatch object
    """

    times = [image[TIME_KEY] for image in image_batch]
    min_time = min(times)
    first_image = image_batch[times.index(min_time)]

    # In case one of the dithers is mis-named, we'll use the most common name
    names = [x[TARGET_KEY] for x in image_batch]
    most_common_name = max(set(names), key=names.count)

    for image in image_batch:
        image["OBJRA"] = first_image["OBJRA"]
        image["OBJDEC"] = first_image["OBJDEC"]
        image[TARGET_KEY] = most_common_name

    return image_batch
