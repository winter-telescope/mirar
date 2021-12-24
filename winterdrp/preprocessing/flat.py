import astropy.io.fits
import numpy as np
import os
from winterdrp.io import create_fits
import logging
import pandas as pd
from collections.abc import Callable
from winterdrp.preprocessing.base_processor import BaseProcessor
from winterdrp.paths import cal_output_dir
# from winterdrp.pipelines.base_pipeline import Pipeline

logger = logging.getLogger(__name__)


class FlatCalibrator(BaseProcessor):

    base_name = "master_flat"
    base_key = "flat"

    def __init__(
            self,
            instrument_vars: dict,
            *args,
            **kwargs
    ):
        super().__init__(instrument_vars, *args, **kwargs)
        self.x_min = instrument_vars["x_min"]
        self.x_max = instrument_vars["x_max"]
        self.y_min = instrument_vars["y_min"]
        self.y_max = instrument_vars["y_max"]
        self.flat_nan_threshold = instrument_vars["flat_nan_threshold"]

    def get_file_path(
            self,
            header: astropy.io.fits.Header,
            sub_dir: str = ""
    ):
        cal_dir = cal_output_dir(sub_dir=sub_dir)
        filtername = header['FILTER']
        name = f"{self.base_name}_{filtername}.fits"
        return os.path.join(cal_dir, name)

    def _apply_to_images(
            self,
            images: list,
            sub_dir: str = ""
    ):
        for i, img in enumerate(images):
            data = img[0].data
            header = img[0].header
            master_flat = self.load_cache_file(self.get_file_path(header, sub_dir=sub_dir))
            if np.any(master_flat < self.flat_nan_threshold):
                master_flat[master_flat < self.flat_nan_threshold] = np.nan
            img[0].data = data / master_flat[0].data
            images[i] = img
        return images

    def make_cache_files(
            self,
            image_list: list,
            preceding_steps: list,
            sub_dir: str = "",
            *args,
            **kwargs
    ):

        logger.info(f'Found {len(image_list)} flat frames')

        with self.open_fits(image_list[0]) as img:
            header = img[0].header

        nx = header['NAXIS1']
        ny = header['NAXIS2']

        filter_list = []

        for flat in image_list:
            with self.open_fits(flat) as img:
                header = img[0].header
            filter_list.append(header['FILTER'])

        image_list = np.array(image_list)

        for filter in list(set(filter_list)):

            mask = np.array([x == filter for x in image_list])

            cut_flat_list = image_list[mask]

            n_frames = np.sum(mask)

            logger.info(f'Found {n_frames} frames for filer {filter}')

            flats = np.zeros((ny, nx, n_frames))

            for i, flat in enumerate(cut_flat_list):
                logger.debug(f'Reading flat {i + 1}/{n_frames}')

                with self.open_fits(flat) as img:

                    # Iteratively apply corrections
                    for f in preceding_steps:
                        img = f(list(img))

                data = np.array([x[0].data for x in list(img)])

                median = np.nanmedian(data[self.x_min:self.x_max, self.y_min:self.y_max])
                flats[:, :, i] = data / median

            logger.info(f'Median combining {n_frames} flats')

            master_flat = np.nanmedian(flats, axis=2)

            with self.open_fits(image_list[0]) as img:
                primary_header = img[0].header

            proc_hdu = create_fits(master_flat, header=primary_header, history='Stacked flat-fielded')
            # Create a new HDU with the processed image data

            primary_header['BZERO'] = 0

            master_flat_path = self.get_file_path(header, sub_dir=sub_dir)

            logger.info(f"Saving stacked 'master flat' for filter {filter} to {master_flat_path}")

            self.save_fits(proc_hdu, master_flat_path)


class StandardFlatCalibrator(FlatCalibrator):

    def __init__(
            self,
            open_fits: Callable[[str], astropy.io.fits.HDUList],
            x_min: float = 0.,
            x_max: float = np.inf,
            y_min: float = 0.,
            y_max: float = np.inf,
            flat_nan_threshold: float = np.nan,
            standard_flat_dir: str = None
    ):
        FlatCalibrator.__init__(
            self,
            open_fits=open_fits,
            x_min=x_min,
            x_max=x_max,
            y_min=y_min,
            y_max=y_max,
            flat_nan_threshold=flat_nan_threshold
        )
        if standard_flat_dir is None:
            err = "For StandardFlatCalibrator, you must specify "
            logger.error(err)
            raise ValueError(err)

    def make_cache_files(
            self,
            image_list: list,
            sub_dir: str = "",
            subtract_bias: Callable[[astropy.io.fits.HDUList], astropy.io.fits.HDUList] = None,
            subtract_dark: Callable[[astropy.io.fits.HDUList], astropy.io.fits.HDUList] = None,
            **kwargs
    ):
        pass

    def get_file_path(
            self,
            header: astropy.io.fits.Header,
            sub_dir: str = ""
    ):
        cal_dir = cal_output_dir(sub_dir=sub_dir)
        filtername = header['FILTER']
        name = f"{self.base_name}_{filtername}.fits"
        return os.path.join(cal_dir, name)