"""
Module containing processors for flat calibration
"""
import logging
import sys
from collections.abc import Callable

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.paths import FLAT_FRAME_KEY, LATEST_SAVE_KEY
from mirar.processors.base_processor import ProcessorPremadeCache, ProcessorWithCache
from mirar.processors.utils.image_selector import select_from_images

logger = logging.getLogger(__name__)


def default_select_flat(
    images: ImageBatch,
) -> ImageBatch:
    """
    Select images tagged as flat

    :param images: set of images
    :return: subset of flat images
    """
    return select_from_images(images, target_values="flat")


class FlatCalibrator(ProcessorWithCache):
    """
    Processor to apply flat calibration
    """

    base_key = "flat"

    def __init__(
        self,
        *args,
        x_min: int = 0,
        x_max: int = sys.maxsize,
        y_min: int = 0,
        y_max: int = sys.maxsize,
        flat_nan_threshold: float = 0.0,
        select_flat_images: Callable[[ImageBatch], ImageBatch] = default_select_flat,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.flat_nan_threshold = flat_nan_threshold
        self.select_cache_images = select_flat_images

    def __str__(self) -> str:
        return "Creates a flat image, divides other images by this image."

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        master_flat = self.get_cache_file(batch)
        master_flat_data = master_flat.get_data()

        mask = master_flat_data <= self.flat_nan_threshold

        if np.sum(mask) > 0:
            master_flat_data[mask] = np.nan

        for image in batch:
            data = image.get_data()
            data = data / master_flat_data
            image.set_data(data)
            image[FLAT_FRAME_KEY] = master_flat[LATEST_SAVE_KEY]

        return batch

    def make_image(
        self,
        images: ImageBatch,
    ) -> Image:
        images = self.select_cache_images(images)

        n_frames = len(images)
        if n_frames == 0:
            err = f"Found {n_frames} suitable flats in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].get_data().shape

        flats = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            median = np.nanmedian(
                img.get_data()[self.x_min : self.x_max, self.y_min : self.y_max]
            )
            flats[:, :, i] = img.get_data() / median

        logger.info(f"Median combining {n_frames} flats")
        master_flat = np.nanmedian(flats, axis=2)

        return Image(master_flat, header=images[0].get_header())


class SkyFlatCalibrator(FlatCalibrator):
    """
    Processor to do flat calibration using sky flats
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, select_flat_images=self.select_sky_flat)

    @staticmethod
    def select_sky_flat(
        images: ImageBatch,
    ) -> ImageBatch:
        """
        Select science images to use as sky flats

        :param images: set of images
        :return: subset of 'sky' images
        """
        return select_from_images(images, key="obsclass", target_values="science")

    def __str__(self) -> str:
        return (
            "Processor to create a sky flat image, divides other images by this image."
        )


class MasterFlatCalibrator(ProcessorPremadeCache, FlatCalibrator):
    """Processor to do flat calibration with a master flate"""
