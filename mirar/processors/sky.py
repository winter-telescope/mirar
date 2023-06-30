"""
Module for sky subtraction
"""
import logging

import numpy as np

from mirar.data import ImageBatch
from mirar.paths import SATURATE_KEY
from mirar.processors.base_processor import ProcessorPremadeCache
from mirar.processors.flat import SkyFlatCalibrator

logger = logging.getLogger(__name__)


class NightSkyMedianCalibrator(SkyFlatCalibrator):
    """
    Processor for sky subtraction
    """

    base_key = "sky"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        master_sky = self.get_cache_file(batch)

        mask = master_sky.get_data() <= self.flat_nan_threshold

        if np.sum(mask) > 0:
            data = master_sky.get_data()
            data[mask] = np.nan
            master_sky.set_data(data)

        for image in batch:
            data = image.get_data()
            header = image.get_header()

            subtract_median = np.nanmedian(data)
            data = data - subtract_median * master_sky.get_data()

            header.append(
                ("SKMEDSUB", subtract_median, "Median sky level subtracted"), end=True
            )
            if SATURATE_KEY in image.header:
                # image[SATURATE_KEY] -= subtract_median
                image[SATURATE_KEY] = 25000

            image.set_data(data)
            image.set_header(header)

        return batch

    def __str__(self) -> str:
        return (
            "Processor to create a median sky background map,"
            "and subtract this from images."
        )


class MasterSkyCalibrator(ProcessorPremadeCache, NightSkyMedianCalibrator):
    """
    Processor to subtract a master sky image
    """
