import numpy as np
import logging
from winterdrp.processors.base_processor import ProcessorWithCache, ProcessorPremadeCache
from collections.abc import Callable
import astropy.io.fits
from winterdrp.processors.utils.image_selector import select_from_images
from winterdrp.paths import latest_save_key, bias_frame_key
from winterdrp.errors import ImageNotFoundError
from winterdrp.data import Image, ImageBatch

logger = logging.getLogger(__name__)


def default_select_bias(
        images: ImageBatch,
) -> ImageBatch:
    return select_from_images(images, target_values="bias")


class BiasCalibrator(ProcessorWithCache):

    base_key = "bias"

    def __init__(
            self,
            select_bias_images: Callable[[ImageBatch], ImageBatch] = default_select_bias,
            *args,
            **kwargs
    ):
        super(BiasCalibrator, self).__init__(*args, **kwargs)
        self.select_cache_images = select_bias_images

    def __str__(self) -> str:
        return f"Creates a bias image, and subtracts this from the other images."

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:

        master_bias = self.get_cache_file(batch)

        for i, image in enumerate(batch):
            data = image.get_data()
            data = data - master_bias.get_data()
            image.set_data(data)
            image[bias_frame_key] = master_bias[latest_save_key]

        return batch

    def make_image(
            self,
            image_batch: ImageBatch,
    ) -> Image:

        images = self.select_cache_images(image_batch)

        n_frames = len(images)

        if n_frames == 0:
            err = f"Found {n_frames} suitable biases in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].get_data().shape

        biases = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            biases[:, :, i] = img.get_data()

        logger.info(f'Median combining {n_frames} biases')
        master_bias = Image(np.nanmedian(biases, axis=2), header=images[0].get_header())

        return master_bias


class MasterBiasCalibrator(ProcessorPremadeCache, BiasCalibrator):

    def __str__(self) -> str:
        return f"Loads a master bias image, and subtracts this from the other images."
