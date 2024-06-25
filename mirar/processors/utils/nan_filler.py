"""
Module for filling NaN values in an image
"""

import logging

import numpy as np

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class NanFiller(BaseImageProcessor):
    """
    Processor to fill image nans
    """

    base_key = "nanfill"

    def description(self):
        return "Processor to pad masked image regions with the median."

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        for image in batch:
            data = image.get_data()
            data[np.isnan(data)] = np.nanmedian(data)
            image.set_data(data)

        return batch
