"""
Module for applying dark corrections
"""
import logging
from collections.abc import Callable

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    EXPTIME_KEY,
    OBSCLASS_KEY,
    SATURATE_KEY,
    STACKED_COMPONENT_IMAGES_KEY,
)
from mirar.processors.base_processor import ProcessorPremadeCache, ProcessorWithCache
from mirar.processors.utils.image_selector import select_from_images

logger = logging.getLogger(__name__)


class MissingDarkError(ImageNotFoundError):
    """
    Error for when a dark image is missing
    """


def default_select_dark(
    images: ImageBatch,
) -> ImageBatch:
    """
    Function to select images in a batch tagged as 'dark'

    :param images: images to filter
    :return: batch of dark images
    """
    return select_from_images(images, key=OBSCLASS_KEY, target_values="dark")


class DarkCalibrator(ProcessorWithCache):
    """
    Processor for applying dark correction
    """

    base_name = "master_dark"
    base_key = "dark"

    def __init__(
        self,
        *args,
        select_cache_images: Callable[[ImageBatch], ImageBatch] = default_select_dark,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.select_cache_images = select_cache_images

    def __str__(self) -> str:
        return (
            "Processor to create a dark image, "
            "and subtracts this from the other images."
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        master_dark = self.get_cache_file(batch)

        for image in batch:
            data = image.get_data()
            data = data - (master_dark.get_data() * image[EXPTIME_KEY])
            image.set_data(data)

            if SATURATE_KEY in image.header:
                image[SATURATE_KEY] -= (
                    np.nanmedian(master_dark.get_data()) * image[EXPTIME_KEY]
                )
        return batch

    def make_image(
        self,
        images: ImageBatch,
    ) -> Image:
        dark_images = self.select_cache_images(images)

        n_frames = len(dark_images)
        if n_frames == 0:
            err = (
                f"Found {n_frames} suitable darks in batch. "
                f"First image in batch has exposure time {images[0][EXPTIME_KEY]}"
                f"and name {images[0][BASE_NAME_KEY]}"
            )
            logger.error(err)
            raise MissingDarkError(err)

        nx, ny = dark_images[0].get_data().shape

        darks = np.zeros((nx, ny, n_frames))

        individual_dark_exptimes, imagenames_key = [], []
        for i, img in enumerate(dark_images):
            dark_exptime = img[EXPTIME_KEY]
            darks[:, :, i] = img.get_data() / dark_exptime
            individual_dark_exptimes.append(str(dark_exptime))
            imagenames_key.append(img[BASE_NAME_KEY])

        logger.debug(f"Median combining {n_frames} darks")
        master_dark_header = dark_images[0].get_header()
        master_dark_header[EXPTIME_KEY] = 1.0
        master_dark_header[COADD_KEY] = n_frames
        master_dark_header["INDIVEXP"] = ",".join(individual_dark_exptimes)
        master_dark_header[STACKED_COMPONENT_IMAGES_KEY] = ",".join(imagenames_key)
        master_dark = Image(np.nanmedian(darks, axis=2), header=master_dark_header)

        return master_dark


class MasterDarkCalibrator(ProcessorPremadeCache, DarkCalibrator):
    """
    Processor to apply master dark corrections
    """
