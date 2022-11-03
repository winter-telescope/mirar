import copy
import astropy.io.fits
import numpy as np
import os
import logging
import pandas as pd
from collections.abc import Callable
from winterdrp.processors.base_processor import ProcessorWithCache, ProcessorPremadeCache
from winterdrp.paths import base_name_key
from winterdrp.processors.utils.image_selector import select_from_images
from winterdrp.errors import ImageNotFoundError


logger = logging.getLogger(__name__)


def default_select_flat(
        images: list[np.ndarray],
        headers: list[astropy.io.fits.Header],
) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
    return select_from_images(images, headers, target_values="dark")


class DarkCalibrator(ProcessorWithCache):
    base_name = "master_dark"
    base_key = "dark"

    def __init__(
            self,
            select_cache_images: Callable[[list, list], [list, list]] = default_select_flat,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.select_cache_images = select_cache_images

    def __str__(self) -> str:
        return f"Processor to create a dark image, and subtracts this from the other images."

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        master_dark, _ = self.get_cache_file(images, headers)

        for i, data in enumerate(images):
            header = headers[i]
            data = data - (master_dark * header["EXPTIME"])
            images[i] = data
            headers[i] = header

        return images, headers

    def make_image(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ):
        images, headers = self.select_cache_images(images, headers)

        n_frames = len(images)
        if n_frames == 0:
            err = f"Found {n_frames} suitable darks in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].shape

        darks = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            dark_exptime = headers[i]['EXPTIME']
            darks[:, :, i] = img / dark_exptime

        logger.info(f'Median combining {n_frames} darks')
        master_dark = np.nanmedian(darks, axis=2)

        return master_dark, headers[0]


class MasterDarkCalibrator(ProcessorPremadeCache, DarkCalibrator):
    pass
