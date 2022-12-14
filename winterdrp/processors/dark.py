import copy
import logging
import os
from collections.abc import Callable

import astropy.io.fits
import numpy as np
import pandas as pd

from winterdrp.data import Image, ImageBatch
from winterdrp.errors import ImageNotFoundError
from winterdrp.paths import BASE_NAME_KEY
from winterdrp.processors.base_processor import (
    ProcessorPremadeCache,
    ProcessorWithCache,
)
from winterdrp.processors.utils.image_selector import select_from_images

logger = logging.getLogger(__name__)


def default_select_flat(
    images: ImageBatch,
) -> ImageBatch:
    return select_from_images(images, target_values="dark")


class DarkCalibrator(ProcessorWithCache):
    base_name = "master_dark"
    base_key = "dark"

    def __init__(
        self,
        select_cache_images: Callable[[ImageBatch], ImageBatch] = default_select_flat,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.select_cache_images = select_cache_images

    def __str__(self) -> str:
        return f"Processor to create a dark image, and subtracts this from the other images."

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        master_dark = self.get_cache_file(batch)

        for i, image in enumerate(batch):
            data = image.get_data()
            data = data - (master_dark.get_data() * image["EXPTIME"])
            image.set_data(data)

        return batch

    def make_image(
        self,
        images: ImageBatch,
    ) -> Image:
        images = self.select_cache_images(images)

        n_frames = len(images)
        if n_frames == 0:
            err = f"Found {n_frames} suitable darks in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].get_data().shape

        darks = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            dark_exptime = img["EXPTIME"]
            darks[:, :, i] = img.get_data() / dark_exptime

        logger.info(f"Median combining {n_frames} darks")
        master_dark = Image(np.nanmedian(darks, axis=2), header=images[0].get_header())

        return master_dark


class MasterDarkCalibrator(ProcessorPremadeCache, DarkCalibrator):
    pass
