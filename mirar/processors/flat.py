"""
Module containing processors for flat calibration
"""
import logging
import os.path
import sys
from collections.abc import Callable

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.paths import BASE_NAME_KEY, FLAT_FRAME_KEY, LATEST_SAVE_KEY
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
        flat_mask_key: str = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.flat_nan_threshold = flat_nan_threshold
        self.select_cache_images = select_flat_images
        self.flat_mask_key = flat_mask_key

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

        logger.debug(f"Found {len(images)} suitable flats in batch")

        n_frames = len(images)
        if n_frames == 0:
            err = f"Found {n_frames} suitable flats in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].get_data().shape

        flats = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            data = img.get_data()

            if self.flat_mask_key is not None:
                if self.flat_mask_key not in img.header.keys():
                    err = (
                        f"Image {img} does not have a mask with key "
                        f"{self.flat_mask_key}"
                    )
                    logger.error(err)
                    raise KeyError(err)

                mask_file = img[self.flat_mask_key]
                logger.debug(f"Masking flat {img[BASE_NAME_KEY]} with mask {mask_file}")
                if not os.path.exists(mask_file):
                    err = f"Mask file {mask_file} does not exist"
                    logger.error(err)
                    raise FileNotFoundError(err)

                mask_img = self.open_fits(mask_file)
                pixels_to_keep = mask_img.get_data().astype(bool)
                mask = ~pixels_to_keep
                logger.debug(
                    f"Masking {np.sum(mask)} pixels in flat {img[BASE_NAME_KEY]}"
                )
                data[mask] = np.nan

            median = np.nanmedian(
                data[self.x_min : self.x_max, self.y_min : self.y_max]
            )
            flats[:, :, i] = data / median

        logger.debug(f"Median combining {n_frames} flats")

        master_flat = np.nanmedian(flats, axis=2)

        return Image(master_flat, header=images[0].get_header())


class SkyFlatCalibrator(FlatCalibrator):
    """
    Processor to do flat calibration using sky flats
    """

    def __init__(self, flat_mask_key=None, *args, **kwargs):
        super().__init__(
            *args,
            **kwargs,
            select_flat_images=self.select_sky_flat,
            flat_mask_key=flat_mask_key,
        )

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
