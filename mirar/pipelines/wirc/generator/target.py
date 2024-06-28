"""
Module to annotate target coordinates on images.
"""

from astropy.coordinates import Angle

from mirar.data import ImageBatch
from mirar.paths import TIME_KEY


def annotate_target_coordinates(image_batch: ImageBatch) -> ImageBatch:
    """
    Function to annotate target coordinates on images.
    For WIRC, this should be the value of RA/DEC in the header of the first image.

    :param image_batch: ImageBatch object
    :return: ImageBatch object
    """

    times = [image[TIME_KEY] for image in image_batch]
    min_time = min(times)
    first_image = image_batch[times.index(min_time)]

    target_ra = Angle(first_image["RA"], unit="hourangle").degree
    target_dec = Angle(first_image["DEC"], unit="degree").degree

    for image in image_batch:
        image["TARGRA"] = target_ra
        image["TARGDEC"] = target_dec

    return image_batch
