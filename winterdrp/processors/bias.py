import numpy as np
import logging
from winterdrp.processors.base_processor import ProcessorWithCache, ProcessorPremadeCache
from collections.abc import Callable
import astropy.io.fits
from winterdrp.processors.utils.image_selector import select_from_images
from winterdrp.paths import latest_save_key, bias_frame_key
from winterdrp.errors import ImageNotFoundError

logger = logging.getLogger(__name__)


def default_select_bias(
        images: list[np.ndarray],
        headers: list[astropy.io.fits.Header],
) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
    return select_from_images(images, headers, target_values="bias")


class BiasCalibrator(ProcessorWithCache):

    base_key = "bias"

    def __init__(
            self,
            select_bias_images: Callable[[list, list], [list, list]] = default_select_bias,
            *args,
            **kwargs
    ):
        super(BiasCalibrator, self).__init__(*args, **kwargs)
        self.select_cache_images = select_bias_images

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        master_bias, master_bias_header = self.get_cache_file(images, headers)

        for i, data in enumerate(images):
            data = data - master_bias
            images[i] = data
            headers[i][bias_frame_key] = master_bias_header[latest_save_key]

        return images, headers

    def make_image(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ):
        images, headers = self.select_cache_images(images, headers)

        n_frames = len(images)
        if n_frames == 0:
            err = f"Found {n_frames} suitable biases in batch"
            logger.error(err)
            raise ImageNotFoundError(err)

        nx, ny = images[0].shape

        biases = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            biases[:, :, i] = img

        logger.info(f'Median combining {n_frames} biases')
        master_bias = np.nanmedian(biases, axis=2)

        return master_bias, headers[0]


class MasterBiasCalibrator(ProcessorPremadeCache, BiasCalibrator):
    pass
