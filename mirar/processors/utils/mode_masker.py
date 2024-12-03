"""
Module for filling NaN values into an image
"""

import logging

import numpy as np
from scipy.stats import mode

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class ModeMasker(BaseImageProcessor):
    """
    Processor to mask the most common value in an image
    """

    base_key = "modemask"

    def description(self):
        return "Processor to mask image pixels with very common values"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        for image in batch:
            data = image.get_data()

            n_pixels = data.size

            mode_value = mode(data[~np.isnan(data)])[0]
            frac = np.sum(data == mode_value) / n_pixels

            while frac > 0.0001:
                data[data == mode_value] = np.nan
                mode_value = mode(data[~np.isnan(data)])[0]
                frac = np.sum(data == mode_value) / n_pixels

            image.set_data(data)

        return batch
