"""
Module containing processors for bias correction
"""
import logging
from collections.abc import Callable

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.paths import BIAS_FRAME_KEY, LATEST_SAVE_KEY, SATURATE_KEY
from mirar.processors.base_processor import ProcessorPremadeCache, ProcessorWithCache
from mirar.processors.utils.image_selector import select_from_images

logger = logging.getLogger(__name__)


def default_select_bias(
    images: ImageBatch,
) -> ImageBatch:
    """
    Returns images in a batch with are tagged as bias

    :param images: set of images
    :return: subset of bias images
    """
    return select_from_images(images, target_values="bias")


class BiasCalibrator(ProcessorWithCache):
    """
    Processor to bias-correct images
    """

    base_key = "bias"

    def __init__(
        self,
        *args,
        select_bias_images: Callable[[ImageBatch], ImageBatch] = default_select_bias,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.select_cache_images = select_bias_images

    def __str__(self) -> str:
        return "Creates a bias image, and subtracts this from the other images."

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        master_bias = self.get_cache_file(batch)

        for image in batch:
            data = image.get_data()
            data = data - master_bias.get_data()
            image.set_data(data)
            image[BIAS_FRAME_KEY] = master_bias[LATEST_SAVE_KEY]
            if SATURATE_KEY in image.header:
                image[SATURATE_KEY] -= np.nanmedian(master_bias.get_data())
        return batch

    def make_image(
        self,
        images: ImageBatch,
    ) -> Image:
        images = self.select_cache_images(images)

        n_frames = len(images)

        if n_frames == 0:
            err = f"Found {n_frames} suitable biases in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].get_data().shape

        biases = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            biases[:, :, i] = img.get_data()

        logger.debug(f"Median combining {n_frames} biases")
        master_bias = Image(np.nanmedian(biases, axis=2), header=images[0].get_header())

        return master_bias


class MasterBiasCalibrator(ProcessorPremadeCache, BiasCalibrator):
    """
    Processor to bias-correct with a master bias
    """

    def __str__(self) -> str:
        return "Loads a master bias image, and subtracts this from the other images."
