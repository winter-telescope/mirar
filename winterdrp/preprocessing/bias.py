from winterdrp.io import create_fits
import os
import numpy as np
import logging
from winterdrp.preprocessing.base_processor import BaseProcessor
from winterdrp.paths import cal_output_dir


logger = logging.getLogger(__name__)


class BiasCalibrator(BaseProcessor):

    base_name = 'master_bias.fits'
    base_key = "bias"

    def get_file_path(self, header, sub_dir=""):
        cal_dir = cal_output_dir(sub_dir=sub_dir)
        return os.path.join(cal_dir, self.base_name)

    def apply_to_images(self, images, sub_dir=""):
        data = images[0].data
        header = images[0].header
        master_bias = self.load_cache_file(self.get_file_path(header, sub_dir=sub_dir))
        images[0].data = data - master_bias[0].data
        return images

    def make_cache_files(
            self,
            image_list: list,
            preceding_steps: list,
            sub_dir: str = "",
            *args,
            **kwargs
    ):

        image_list = image_list

        logger.info(f'Found {len(image_list)} bias frames')

        with self.open_fits(image_list[0]) as img:
            header = img[0].header

        nx = header['NAXIS1']
        ny = header['NAXIS2']

        nframes = len(image_list)

        biases = np.zeros((ny, nx, nframes))

        for i, bias in enumerate(image_list):
            logger.debug(f'Reading bias {i + 1}/{nframes}')
            with self.open_fits(bias) as img:

                # Iteratively apply corrections
                for f in preceding_steps:
                    img = f(img)

            biases[:, :, i] = np.array([x[0].data for x in img])

        logger.info(f'Median combining {nframes} biases')

        master_bias = np.nanmedian(biases, axis=2)

        with self.open_fits(image_list[0]) as img:
            primary_header = img[0].header

        proc_hdu = create_fits(master_bias, header, history='median stacked bias')
        # Create a new HDU with the processed image data
        proc_hdu.header = primary_header  # Copy over the header from the raw file

        master_bias_path = self.get_file_path(header, sub_dir=sub_dir)
        logger.info(f"Saving stacked 'master bias' to {master_bias_path}")
        self.save_fits(proc_hdu, master_bias_path)
