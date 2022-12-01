import logging

import astropy.io.fits
import numpy as np

from winterdrp.data import ImageBatch
from winterdrp.processors.base_processor import ProcessorPremadeCache
from winterdrp.processors.flat import SkyFlatCalibrator

logger = logging.getLogger(__name__)


class NightSkyMedianCalibrator(SkyFlatCalibrator):

    base_key = "sky"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        master_sky = self.get_cache_file(batch)

        mask = master_sky.get_data() <= self.flat_nan_threshold

        if np.sum(mask) > 0:
            master_sky[mask] = np.nan

        for i, image in enumerate(batch):
            data = image.get_data()
            header = image.get_header()

            subtract_median = np.nanmedian(data)
            data = data - subtract_median * master_sky.get_data()

            header.append(
                ("SKMEDSUB", subtract_median, "Median sky level subtracted"), end=True
            )

            image.set_data(data)
            image.set_header(header)

        return batch

    def __str__(self) -> str:
        return f"Processor to create a median sky background map, and subtract this from images."


class MasterSkyCalibrator(ProcessorPremadeCache, NightSkyMedianCalibrator):
    pass
