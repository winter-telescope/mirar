from winterdrp.io import create_fits
import os
import numpy as np
import logging
from winterdrp.processors.base_processor import ProcessorWithCache
from winterdrp.paths import cal_output_dir
from collections.abc import Callable
import pandas as pd
import astropy.io.fits


logger = logging.getLogger(__name__)


def default_select_bias(
       observing_log: pd.DataFrame
) -> [str]:
    mask = observing_log["OBSTYPE"].lower() == "bias"
    return list(observing_log[mask]["RAWIMAGEPATH"])


class BiasCalibrator(ProcessorWithCache):

    base_name = 'master_bias.fits'
    base_key = "bias"

    def __init__(
            self,
            select_cache_images: Callable[[pd.DataFrame], list] = None,
            *args,
            **kwargs
    ):
        super(BiasCalibrator, self).__init__(*args, **kwargs)
        if select_cache_images is None:

            def select_cache_images(x):
                return self.select_from_log(x, "bias")

        self.select_cache_images = select_cache_images

    def get_file_path(self, sub_dir=""):
        cal_dir = cal_output_dir(sub_dir=sub_dir)
        return os.path.join(cal_dir, self.base_name)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, data in enumerate(images):
            header = headers[i]
            master_bias, _ = self.load_cache_file(self.get_file_path(sub_dir=self.night_sub_dir))
            data = data - master_bias
            # header["CALSTEPS"] += "bias,"
            images[i] = data
            headers[i] = header

        return images, headers

    def make_cache_files(
            self,
            image_list: list[str],
    ):

        image_list = image_list

        logger.info(f'Found {len(image_list)} bias frames')

        _, primary_header = self.open_fits(image_list[0])

        nx = primary_header['NAXIS1']
        ny = primary_header['NAXIS2']

        nframes = len(image_list)

        biases = np.zeros((ny, nx, nframes))

        for i, bias in enumerate(image_list):
            logger.debug(f'Reading bias {i + 1}/{nframes}')
            img, header = self.open_fits(bias)

            # Iteratively apply corrections
            # for f in preceding_steps:
            #    img, header = f(list(img), list(header), sub_dir=sub_dir)

            biases[:, :, i] = img

        logger.info(f'Median combining {nframes} biases')

        master_bias = np.nanmedian(biases, axis=2)

        proc_hdu = create_fits(master_bias, primary_header)
        # Create a new HDU with the processed image data
        proc_hdu.header = primary_header  # Copy over the header from the raw file

        master_bias_path = self.get_file_path(sub_dir=self.night_sub_dir)
        try:
            os.makedirs(os.path.dirname(master_bias_path))
        except OSError:
            pass

        logger.info(f"Saving stacked 'master bias' to {master_bias_path}")
        self.save_fits(master_bias, primary_header, master_bias_path)
