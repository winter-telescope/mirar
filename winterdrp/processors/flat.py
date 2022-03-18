import astropy.io.fits
import numpy as np
import os
import logging
import pandas as pd
import sys
from collections.abc import Callable
import copy
from winterdrp.processors.base_processor import BaseProcessor, ProcessorWithCache
from winterdrp.paths import cal_output_dir

logger = logging.getLogger(__name__)


class FlatCalibrator(ProcessorWithCache):

    base_name = "master_flat"
    base_key = "flat"

    def __init__(
            self,
            x_min: int = 0,
            x_max: int = sys.maxsize,
            y_min: int = 0,
            y_max: int = sys.maxsize,
            flat_nan_threshold: float = np.nan,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.flat_nan_threshold = flat_nan_threshold

    def select_cache_images(self, x):
        return self.select_from_log(x, "flat")

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, data in enumerate(images):
            header = headers[i]
            flat, _ = self.load_cache_file(self.get_file_path(header))
            if np.any(~(flat > self.flat_nan_threshold)):
                flat[~(flat > self.flat_nan_threshold)] = np.nan

            data = data / flat
            header["CALSTEPS"] += "flat,"
            images[i] = data
            headers[i] = header
        return images, headers

    def get_file_path(
            self,
            header: astropy.io.fits.Header,
    ) -> str:
        cal_dir = cal_output_dir(sub_dir=self.night_sub_dir)

        filtername = header['FILTER'].replace(" ", "_")
        name = f"{self.base_name}_{filtername}.fits"
        return os.path.join(cal_dir, name)

    def make_cache_files(
            self,
            image_paths: list[str],
    ):

        logger.info(f'Found {len(image_paths)} flat frames')

        _, primary_header = self.open_fits(image_paths[0])

        nx = primary_header['NAXIS1']
        ny = primary_header['NAXIS2']

        filter_list = []

        for flat in image_paths:
            _, header = self.open_fits(flat)
            filter_list.append(header['FILTER'])

        image_paths = np.array(image_paths)

        for filt in list(set(filter_list)):

            mask = np.array([x == filt for x in filter_list])

            cut_flat_list = image_paths[mask]

            n_frames = np.sum(mask)

            logger.info(f'Found {n_frames} frames for filer {filt}')

            flats = np.zeros((ny, nx, n_frames))

            for i, flat in enumerate(cut_flat_list):
                logger.debug(f'Reading flat {i + 1}/{n_frames}')

                img, header = self.open_fits(flat)

                # Iteratively apply corrections
                for p in self.preceding_steps:
                    img, header = p.apply([img], [header])

                median = np.nanmedian(img[0][self.x_min:self.x_max, self.y_min:self.y_max])
                flats[:, :, i] = img / median

            logger.info(f'Median combining {n_frames} flats')

            master_flat = np.nanmedian(flats, axis=2)

            # Create a new HDU with the processed image data

            primary_header['BZERO'] = 0
            primary_header["FILTER"] = filt
            primary_header['OBJECT'] = "flat"

            master_flat_path = self.get_file_path(primary_header)

            logger.info(f"Saving stacked 'master flat' for f {filt} to {master_flat_path}")

            self.save_fits(master_flat, primary_header, master_flat_path)


class SkyFlatCalibrator(FlatCalibrator):

    def __init__(
            self,
            use_full_night: bool = True,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)

        self.use_full_night = use_full_night

    def select_cache_images(
            self,
            observing_log: pd.DataFrame
    ) -> list[str]:

        mask = np.logical_and(
            observing_log["NIGHT"] == self.night,
            observing_log["OBSCLASS"] == "science"
        )
        obs = observing_log[mask]

        if self.use_full_night:
            return list(obs["RAWIMAGEPATH"])

        else:
            raise NotImplementedError






