"""
Module containing processors which mask pixels
"""
import logging
from pathlib import Path

import numpy as np

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

# MASK_VALUE = -99.
MASK_VALUE = np.nan


class MaskPixels(BaseImageProcessor):
    """
    Processor to apply bias calibration
    """

    base_key = "mask"

    def __init__(self, mask_path: str | Path):
        super().__init__()
        self.mask = None
        self.mask_path = Path(mask_path)

    def __str__(self) -> str:
        return f"Processor to mask bad pixels using a pre-defined map: {self.mask_path}"

    def get_mask(self):
        """
        loads mask if needed, and returns it

        :return: mask
        """
        if self.mask is None:
            self.mask = self.open_fits(self.mask_path)
        return self.mask

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            data = image.get_data()
            mask = self.get_mask().get_data()
            mask = mask != 0
            data[mask] = MASK_VALUE
            image.set_data(data)

        return batch
