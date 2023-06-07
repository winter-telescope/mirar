"""
Module containing functions and processors to filter images
"""
import logging

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor, CleanupProcessor

logger = logging.getLogger(__name__)


def filter_images(
    images: ImageBatch,
    header_key: str = "target",
    reject_values: str | list[str] = "science",
) -> ImageBatch:
    """
    Finds the subset of images in the batch with do not have a value of <header_key>
    equal to a value in <reject values>

    :param images: images to filter
    :param header_key: key to check
    :param reject_values: unacceptable value(s)
    :return: subset of passing images
    """

    # Enforce string in list for later matching
    if not isinstance(reject_values, list):
        reject_values = [str(reject_values)]
    else:
        reject_values = [str(x) for x in reject_values]

    passing_images = ImageBatch()

    for image in images:
        if str(image[header_key]) not in reject_values:
            passing_images.append(image)

    return passing_images


class ImageRejector(BaseImageProcessor, CleanupProcessor):
    """
    Processor to reject images based on the headers
    """

    base_key = "reject"

    def __init__(self, *args: tuple[str, str | list[str]]):
        super().__init__()
        self.rejects = args

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for header_key, reject_values in self.rejects:
            batch = filter_images(
                batch, header_key=header_key, reject_values=reject_values
            )

        return batch
