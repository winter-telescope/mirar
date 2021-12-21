import logging
import os
import numpy as np
from astropy.io import fits
from winterdrp.preprocessing import BiasCalibrator, DarkCalibrator, FlatCalibrator
from winterdrp.calibrate.sourceextractor import run_sextractor
from winterdrp.io import create_fits
from glob import glob
import pandas as pd
from winterdrp.paths import \
    cal_output_dir,\
    parse_image_list, \
    reduced_img_dir, \
    reduced_img_path, \
    raw_img_dir, \
    astrometry_output_dir, \
    observing_log_dir

logger = logging.getLogger(__name__)


class Pipeline:

    subclasses = {}

    # Set up elements to use
    bias = True
    dark = True
    flat = True
    stack = True
    dither = True
    astrometry = ("GAIA", 9., 13.)
    photometry_cal = dict()

    # Fix keys from header to save in log
    # The first key should sort images in order

    header_keys = list()

    # Pixel range for flat-fielding
    x_min = 0.
    x_max = np.inf
    y_min = 0.
    y_max = np.inf

    def __init__(self, *args, **kwargs):
        self.bias_calibrator = [None, BiasCalibrator(self.open_fits, *args, **kwargs)][self.bias]
        self.dark_calibrator = [None, DarkCalibrator(self.open_fits, *args, **kwargs)][self.dark]
        self.flats_calibrator = [None, FlatCalibrator(
            self.open_fits,
            x_min=self.x_min,
            x_max=self.x_max,
            y_min=self.y_min,
            y_max=self.y_max,
            *args, **kwargs
        )][self.flat]
        self.observing_logs_cache = dict()

    def open_fits(self, path):
        img = fits.open(path)
        img = self.reformat_raw_data(img)
        return img

    @staticmethod
    def reformat_raw_data(img):
        return img

    def make_calibration_files(self, sub_dir):

        cal_dict = parse_image_list(sub_dir)

        cal_dir = cal_output_dir(sub_dir)

        # Make calibration directory, unless it already exists

        try:
            os.makedirs(cal_dir)
        except OSError:
            pass

        logger.info(f"Making calibration files for directory {raw_img_dir(sub_dir)}")

        if self.bias:
            self.bias_calibrator.make_calibration_files(
                cal_dict["bias"],
                cal_dir=cal_dir,
                open_fits=self.open_fits
            )

        if self.dark:
            self.dark_calibrator.make_calibration_files(
                image_list=cal_dict["dark"],
                sub_dir=sub_dir,
                subtract_bias=self.subtract_bias
            )

        if self.flat:
            self.flats_calibrator.make_calibration_files(
                cal_dict["flats"],
                cal_dir=cal_dir,
                open_fits=self.open_fits,
                subtract_bias=self.subtract_bias,
                subtract_dark=self.subtract_dark
            )

    def reduce_image(self, img, sub_dir=""):

        img = self.subtract_bias(img, sub_dir=sub_dir)
        img = self.subtract_dark(img, sub_dir=sub_dir)
        img = self.divide_flat(img, sub_dir=sub_dir)

        return img

    def preprocess_images(self, sub_dir="", raw_image_list=None, reprocess=True):

        if raw_image_list is None:
            raw_image_list = parse_image_list(sub_dir, group_by_object=False)

        # Try making output directory, unless it exists

        output_dir = reduced_img_dir(sub_dir)

        try:
            os.makedirs(output_dir)
        except OSError:
            pass

        nframes = len(raw_image_list)

        proccessed_list = []

        # Loop over science images

        for i, raw_img_path in enumerate(raw_image_list):

            img_name = os.path.basename(raw_img_path)

            logger.debug(f"Processing image {i + 1}/{nframes} ({img_name})")

            output_path = reduced_img_path(img_name, sub_dir=sub_dir)

            if np.logical_and(os.path.exists(output_path), reprocess is False):
                logger.debug(f"Skipping image {img_name}, because it has already "
                             f"been processed and 'reprocess' is False.")
                continue

            with self.open_fits(raw_img_path) as img:
                header = img[0].header

                if header['OBSTYPE'] not in ['science', "object"]:
                    logger.debug(f'Obstype is not science, skipping {raw_img_path}')
                    continue

                data_redux = self.reduce_image(img, sub_dir=sub_dir)[0].data

                proc_hdu = create_fits(data_redux, header=header, history=None)

                proc_hdu.header['BZERO'] = 0

                # Write the reduced frame to disk

                logger.debug(f"Saving processed image to {output_path}")
                proccessed_list.append(output_path)
                proc_hdu.writeto(output_path, overwrite=True)

        return proccessed_list

    def subtract_bias(self, img, sub_dir=""):

        if self.bias:
            img = self.bias_calibrator.apply_calibration(img, sub_dir=sub_dir)

        return img

    def subtract_dark(self, img, sub_dir=""):

        if self.dark:
            img = self.dark_calibrator.apply_calibration(img, sub_dir=sub_dir)

        return img

    def divide_flat(self, img, sub_dir=""):

        if self.flat:
            img = self.flats_calibrator.apply_calibration(img, sub_dir=sub_dir)

        return img

    @staticmethod
    def apply_astrometry(sub_dir="", redux_image_list=None, reprocess=True):

        if redux_image_list is None:
            redux_image_list = parse_image_list(sub_dir, group_by_object=False, base_dir_f=reduced_img_dir)

        # Try making output directory, unless it exists

        output_dir = astrometry_output_dir(sub_dir)

        try:
            os.makedirs(output_dir)
        except OSError:
            pass

        # First run Sextractor

        run_sextractor(
            redux_image_list,
            output_dir=output_dir,
            reprocess=reprocess
        )

    def export_observing_log(self, sub_dir=""):

        log = self.load_observing_log(sub_dir=sub_dir)

        path = self.get_observing_log_path(sub_dir=sub_dir)

        logger.debug(f"Saving to {path}")

        log.to_csv(path)

    @staticmethod
    def get_observing_log_path(sub_dir=""):
        base_dir = observing_log_dir(sub_dir=sub_dir)
        return os.path.join(base_dir, "observing_log.csv")

    def load_observing_log(self, sub_dir=""):

        path = self.get_observing_log_path(sub_dir=sub_dir)

        if path in self.observing_logs_cache:
            log = self.observing_logs_cache[path]
        else:
            log = self.parse_observing_log(sub_dir=sub_dir)
            self.observing_logs_cache[path] = log
            self.export_observing_log(sub_dir=sub_dir)

        return log

    def parse_observing_log(self, sub_dir=""):

        img_list = glob(f'{raw_img_dir(sub_dir)}/*.fits')

        log_data = []

        for img_file in img_list:
            header = self.open_fits(img_file)[0].header

            row = []

            for key in self.header_keys:
                row.append(header[key])

            row.append(img_file)

            log_data.append(row)

        log = pd.DataFrame(log_data, columns=self.header_keys + ["RAWIMAGEPATH"])
        log = log.sort_values(by=self.header_keys[0]).reset_index(drop=True)

        return log



