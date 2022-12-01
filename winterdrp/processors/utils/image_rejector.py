import logging

import astropy.io.fits
import numpy as np

from winterdrp.data import Dataset, ImageBatch
from winterdrp.processors.base_processor import BaseImageProcessor, CleanupProcessor

logger = logging.getLogger(__name__)


def filter_images(
    images: ImageBatch,
    header_key: str = "target",
    reject_values: str | list[str] = "science",
) -> ImageBatch:

    # Enforce string in list for later matching
    if not isinstance(reject_values, list):
        reject_values = [str(reject_values)]
    else:
        reject_values = [str(x) for x in reject_values]

    passing_images = ImageBatch()

    for i, image in enumerate(images):
        if str(image[header_key]) not in reject_values:
            passing_images.append(image)

    return passing_images


class ImageRejector(BaseImageProcessor, CleanupProcessor):

    base_key = "reject"

    def __init__(self, *args: tuple[str, str | list[str]], **kwargs):
        super().__init__(*args, **kwargs)
        self.rejects = args

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        for (header_key, reject_values) in self.rejects:

            batch = filter_images(
                batch, header_key=header_key, reject_values=reject_values
            )

        return batch
