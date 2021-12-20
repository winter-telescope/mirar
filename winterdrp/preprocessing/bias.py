from winterdrp.io import create_fits
import os
import numpy as np
import logging
from winterdrp.preprocessing.base_calibrator import BaseCalibrator
from winterdrp.paths import cal_output_dir


logger = logging.getLogger(__name__)


class BiasCalibrator(BaseCalibrator):

    base_name = 'master_bias.fits'

    def get_file_path(self, header, sub_dir=""):

        cal_dir = cal_output_dir(sub_dir=sub_dir)

        return os.path.join(cal_dir, self.base_name)

    def apply_calibration(self, img, sub_dir=""):
        data = img[0].data
        header = img[0].header
        master_bias = self.load_calibrator_file(self.get_file_path(header, sub_dir=sub_dir))
        img[0].data = data - master_bias[0].data
        return img

    def make_calibration_files(self, image_list, sub_dir="", *args, **kwargs):

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
                biases[:, :, i] = img[0].data

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


# base_mbias_name = 'master_bias.fits'
#
#
# def make_master_bias(bias_list, cal_dir, open_fits):
#
#     if len(bias_list) > 0:
#
#         logger.info(f'Found {len(bias_list)} bias frames')
#
#         with open_fits(bias_list[0]) as img:
#             header = img[0].header
#
#         nx = header['NAXIS1']
#         ny = header['NAXIS2']
#
#         nframes = len(bias_list)
#
#         biases = np.zeros((ny, nx, nframes))
#
#         for i, bias in enumerate(bias_list):
#             logger.debug(f'Reading bias {i+1}/{nframes}')
#             with open_fits(bias) as img:
#                 biases[:, :, i] = img[0].data
#
#         logger.info(f'Median combining {nframes} biases')
#
#         master_bias = np.nanmedian(biases, axis=2)
#
#         with open_fits(bias_list[0]) as img:
#             primary_header = img[0].header
#
#         proc_hdu = create_fits(master_bias)  # Create a new HDU with the processed image data
#         proc_hdu.header = primary_header      # Copy over the header from the raw file
#         proc_hdu.header.add_history('median stacked bias')
#
#         mbias_path = os.path.join(cal_dir, base_mbias_name)
#
#         logger.info(f"Saving stacked 'master bias' to {mbias_path}")
#
#         proc_hdu.writeto(mbias_path, overwrite=True)
#
#     else:
#         logger.debug("Skipping master bias creation.")
#
#
# def load_master_bias(cal_dir, open_fits, header=None):
#
#     # Try to load bias image
#
#     try:
#         with open_fits(os.path.join(cal_dir, base_mbias_name)) as img:
#             master_bias = img[0].data
#
#     except FileNotFoundError:
#
#         try:
#
#             nx = header['NAXIS1']
#             ny = header['NAXIS2']
#
#             master_bias = np.zeros((ny, nx))
#
#             logger.warning("No master bias found. No bias correction will be applied.")
#
#         except KeyError:
#
#             err = "No master bias files found, and no header info provided to create a dummy image."
#             logger.error(err)
#             raise FileNotFoundError(err)
#
#     return master_bias
