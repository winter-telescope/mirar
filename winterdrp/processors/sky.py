import astropy.io.fits
import numpy as np
import os
import logging
import pandas as pd
from collections.abc import Callable
from astropy.time import Time
from winterdrp.processors.base_processor import BaseProcessor, ProcessorWithCache
from winterdrp.paths import cal_output_dir

logger = logging.getLogger(__name__)


class SkyMedianCalibrator(ProcessorWithCache):

    base_key = "sky"

    def __init__(
            self,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, data in enumerate(images):
            header = headers[i]
            sky, _ = self.select_sky(data, header, images, headers)
            data = data - np.nanmedian(sky)
            header["CALSTEPS"] += "sky,"
            images[i] = data
            headers[i] = header
        return images, headers

    def select_sky(
            self,
            data: np.ndarray,
            header: astropy.io.fits.header,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.header]
    ) -> tuple[np.ndarray, astropy.io.fits.header]:

        ref = Time(header["UTCTIME"], format='isot', scale='utc').jd

        deltas = np.array([
            abs(Time(x["UTCTIME"], format='isot', scale='utc').jd - ref)
            for x in headers
        ])

        closest_delta_t = min(deltas[deltas > 0.])

        match_index = list(deltas).index(closest_delta_t)

        return images[match_index], headers[match_index]