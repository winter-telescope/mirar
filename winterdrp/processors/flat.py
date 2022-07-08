import astropy.io.fits
import numpy as np
import logging
import sys
from winterdrp.processors.base_processor import ProcessorWithCache
from winterdrp.processors.utils.image_selector import select_from_images
from collections.abc import Callable

logger = logging.getLogger(__name__)


def default_select_flat(
        images: list[np.ndarray],
        headers: list[astropy.io.fits.Header],
) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
    return select_from_images(images, headers, target_values="flat")


class FlatCalibrator(ProcessorWithCache):

    base_key = "flat"

    def __init__(
            self,
            x_min: int = 0,
            x_max: int = sys.maxsize,
            y_min: int = 0,
            y_max: int = sys.maxsize,
            flat_nan_threshold: float = 0.,
            select_flat_images: Callable[[list, list], [list, list]] = default_select_flat,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.flat_nan_threshold = flat_nan_threshold
        self.select_cache_images = select_flat_images

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        master_flat, _ = self.get_cache_file(images, headers)

        mask = master_flat <= self.flat_nan_threshold

        if np.sum(mask) > 0:
            master_flat[mask] = np.nan

        for i, data in enumerate(images):
            header = headers[i]

            data = data / master_flat

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

        nx, ny = images[0].shape

        flats = np.zeros((nx, ny, n_frames))

        for i, img in enumerate(images):
            median = np.nanmedian(img[self.x_min:self.x_max, self.y_min:self.y_max])
            flats[:, :, i] = img / median

        logger.info(f'Median combining {n_frames} flats')
        master_flat = np.nanmedian(flats, axis=2)

        return master_flat, headers[0]


class SkyFlatCalibrator(FlatCalibrator):

    def __init__(
            self,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs, select_flat_images=self.select_sky_flat)

    @staticmethod
    def select_sky_flat(
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        return select_from_images(images, headers, header_key="obsclass", target_values="science")


# class OldSkyFlatCalibrator(SkyFlatCalibrator):
#
#     def select_cache_images(
#             self,
#             observing_log: pd.DataFrame
#     ) -> list[str]:
#
#         mask = observing_log["OBSCLASS"] == "science"
#         obs = observing_log[mask]
#         nights = list(sorted(set(obs["NIGHT"])))
#
#         if len(nights) == 1:
#             err = "Only found one night in observing log. " \
#                   "Cannot use a previous night for creating flats. " \
#                   "Try adjusting the 'log_history_nights' parameter of the pipeline, " \
#                   "to read in more history."
#
#             logger.error(err)
#             raise Exception(err)
#
#         previous_night = nights[nights.index(self.night) - 1]
#
#         previous = obs[obs["NIGHT"] == previous_night]
#
#         return list(previous["RAWIMAGEPATH"])






