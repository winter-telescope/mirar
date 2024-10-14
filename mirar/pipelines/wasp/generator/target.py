"""
Module to annotate target coordinates on images.
"""

from astropy.coordinates import Angle

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

    target_ra = Angle(first_image["RA"], unit="hourangle").degree
    target_dec = Angle(first_image["DEC"], unit="degree").degree

    # In case one of the dithers is mis-named, we'll use the most common name
    names = [x[TARGET_KEY] for x in image_batch]
    most_common_name = max(set(names), key=names.count)

    for image in image_batch:
        image["OBJRA"] = target_ra
        image["OBJDEC"] = target_dec
        image["RA"] = first_image["RA"]
        image["DEC"] = first_image["DEC"]
        image[TARGET_KEY] = most_common_name

    return image_batch
