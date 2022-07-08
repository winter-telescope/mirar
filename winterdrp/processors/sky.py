import astropy.io.fits
import numpy as np
import os
import logging
import pandas as pd
from collections.abc import Callable
from astropy.time import Time
from winterdrp.processors.base_processor import ProcessorWithCache
from winterdrp.processors.flat import SkyFlatCalibrator
from winterdrp.paths import cal_output_dir, saturate_key

logger = logging.getLogger(__name__)


class NightSkyMedianCalibrator(SkyFlatCalibrator):

    base_key = "sky"

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        master_sky, _ = self.get_cache_file(images, headers)

        mask = master_sky <= self.flat_nan_threshold

        if np.sum(mask) > 0:
            master_sky[mask] = np.nan

        for i, data in enumerate(images):
            header = headers[i]

            subtract_median = np.nanmedian(data)
            data = data - subtract_median * master_sky

            header.append(('SKMEDSUB', subtract_median, 'Median sky level subtracted'), end=True)

            images[i] = data
            headers[i] = header

        return images, headers


# class OldNightSkyMedianCalibrator(
#     NightSkyMedianCalibrator,
#     OldSkyFlatCalibrator
# ):
#     pass

    # def select_sky(
    #         self,
    #         data: np.ndarray,
    #         header: astropy.io.fits.header,
    #         images: list[np.ndarray],
    #         headers: list[astropy.io.fits.header]
    # ) -> tuple[np.ndarray, astropy.io.fits.header]:
    #
    #     ref = Time(header["UTCTIME"], format='isot', scale='utc').jd
    #
    #     deltas = np.array([
    #         abs(Time(x["UTCTIME"], format='isot', scale='utc').jd - ref)
    #         for x in headers
    #     ])
    #
    #     closest_delta_t = min(deltas[deltas > 0.])
    #
    #     match_index = list(deltas).index(closest_delta_t)
    #
    #     return images[match_index], headers[match_index]