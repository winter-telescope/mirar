import logging
import os

import astropy.io.fits
import numpy as np
from astropy.io import fits
from winterdrp.preprocessing import get_processor
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


raw_img_key = "RAWIMAGEPATH"

class Pipeline:

    # Set up elements to use
    astrometry = ("GAIA", 9., 13.)
    photometry_cal = dict()

    image_steps = list()

    # Fix keys from header to save in log
    # The first key should sort images in order

    header_keys = list()

    # Key from log which defines batch grouping
    # By default this is the raw image path, i.e each each is processed separately

    batch_split_keys = ["RAWIMAGEPATH"]

    # Pixel range for flat-fielding
    x_min = 0.
    x_max = np.inf
    y_min = 0.
    y_max = np.inf
    flat_nan_threshold = np.nan
    standard_flats_dir = None

    def __init__(self, *args, **kwargs):
        # self.bias_calibrator = [None, BiasCalibrator(self.open_fits, *args, **kwargs)][self.bias]
        # self.dark_calibrator = [None, DarkCalibrator(self.open_fits, *args, **kwargs)][self.dark]
        # self.flats_calibrator = [None, FlatCalibrator(
        #     self.open_fits,
        #     x_min=self.x_min,
        #     x_max=self.x_max,
        #     y_min=self.y_min,
        #     y_max=self.y_max,
        #     *args, **kwargs
        # )][self.flat]

        instrument_vars = dict([(x, getattr(self, x)) for x in dir(self)])

        self.observing_logs_cache = dict()
        self.processors = dict()

        for name in self.image_steps:
            # if isinstance(processor_args, tuple):
            #     name = processor_args[0]
            #     args += processor_args[1:]
            # else:
            #     name = processor_args
            self.processors[name] = get_processor(name, instrument_vars, *args, **kwargs)

    def open_fits(
            self,
            path: str
    ) -> [np.array, astropy.io.fits.Header]:
        with fits.open(path) as raw:
            img = self.reformat_raw_data(raw)
            data = img[0].data
            header = img[0].header
        return data, header

    @staticmethod
    def reformat_raw_data(img):
        return img

    def make_calibration_files(self, sub_dir=""):

        observing_log = self.load_observing_log(sub_dir=sub_dir)

        cal_dir = cal_output_dir(sub_dir)

        # Make calibration directory, unless it already exists

        try:
            os.makedirs(cal_dir)
        except OSError:
            pass

        logger.info(f"Making calibration files for directory {raw_img_dir(sub_dir)}")

        preceding_steps = []

        for processor_name in self.image_steps:
            processor = self.processors[processor_name]
            image_list = processor.select_cache_images(observing_log)
            if len(image_list) > 0:
                processor.make_cache_files(
                    image_list=image_list,
                    sub_dir=sub_dir,
                    preceding_steps=preceding_steps,
                )
            preceding_steps.append(self.processors[processor_name].apply_to_images)

        # if self.bias:
        #
        #     bias_files = log[log["OBJECT"] == "bias"]["RAWIMAGEPATH"]
        #
        #     self.bias_calibrator.make_cache_files(
        #         image_list=bias_files,
        #         cal_dir=cal_dir,
        #         open_fits=self.open_fits
        #     )
        #
        # if self.dark:
        #
        #     dark_files = log[log["OBJECT"] == "dark"]["RAWIMAGEPATH"]
        #
        #     self.dark_calibrator.make_cache_files(
        #         image_list=dark_files,
        #         sub_dir=sub_dir,
        #         subtract_bias=self.subtract_bias
        #     )
        #
        # if self.flat:
        #
        #     flat_files = log[log["OBJECT"] == "flat"]["RAWIMAGEPATH"]
        #
        #     self.flats_calibrator.make_cache_files(
        #         image_list=flat_files,
        #         cal_dir=cal_dir,
        #         open_fits=self.open_fits,
        #         subtract_bias=self.subtract_bias,
        #         subtract_dark=self.subtract_dark
        #     )

    def split_raw_images_into_batches(
            self,
            sub_dir: str = ""
    ) -> list:
        observing_log = self.load_observing_log(sub_dir=sub_dir)

        mask = observing_log["OBSCLASS"] == "science"
        obs = observing_log[mask]

        split_vals = []

        for row in obs.itertuples():
            sv = ""
            for key in self.batch_split_keys:
                sv += str(getattr(row, key))

            split_vals.append(sv)

        split_vals = np.array(split_vals)

        raw_image_path_batches = [
            list(obs[split_vals == x][raw_img_key])
            for x in list(set(split_vals))
        ]

        return sorted(raw_image_path_batches)

    def reduce_images(
            self,
            images: list = None,
            sub_dir: str = ""
    ) -> list:

        if images is None:
            images = parse_image_list(sub_dir, group_by_object=False)

        for processor in self.image_steps:
            images = self.processors[processor].apply(images, sub_dir=sub_dir)

        return images

    def process_images(self, sub_dir="", raw_image_list=None, reprocess=True):

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
                header = img.header

                if header['OBSTYPE'] not in ['science', "object"]:
                    logger.debug(f'Obstype is not science, skipping {raw_img_path}')
                    continue

                data_redux = np.array(self.reduce_single_image(img, sub_dir=sub_dir))

                proc_hdu = create_fits(data_redux, header=header, history=None)

                proc_hdu.header['BZERO'] = 0

                # Write the reduced frame to disk

                logger.debug(f"Saving processed image to {output_path}")
                proccessed_list.append(output_path)
                proc_hdu.writeto(output_path, overwrite=True)

        return proccessed_list

    # def preprocess_images(self, sub_dir="", raw_image_list=None, reprocess=True):
    #
    #     if raw_image_list is None:
    #         raw_image_list = parse_image_list(sub_dir, group_by_object=False)
    #
    #     # Try making output directory, unless it exists
    #
    #     output_dir = reduced_img_dir(sub_dir)
    #
    #     try:
    #         os.makedirs(output_dir)
    #     except OSError:
    #         pass
    #
    #     nframes = len(raw_image_list)
    #
    #     proccessed_list = []
    #
    #     # Loop over science images
    #
    #     for i, raw_img_path in enumerate(raw_image_list):
    #
    #         img_name = os.path.basename(raw_img_path)
    #
    #         logger.debug(f"Processing image {i + 1}/{nframes} ({img_name})")
    #
    #         output_path = reduced_img_path(img_name, sub_dir=sub_dir)
    #
    #         if np.logical_and(os.path.exists(output_path), reprocess is False):
    #             logger.debug(f"Skipping image {img_name}, because it has already "
    #                          f"been processed and 'reprocess' is False.")
    #             continue
    #
    #         with self.open_fits(raw_img_path) as img:
    #             header = img.header
    #
    #             if header['OBSTYPE'] not in ['science', "object"]:
    #                 logger.debug(f'Obstype is not science, skipping {raw_img_path}')
    #                 continue
    #
    #             data_redux = np.array(self.reduce_single_image(img, sub_dir=sub_dir))
    #
    #             proc_hdu = create_fits(data_redux, header=header, history=None)
    #
    #             proc_hdu.header['BZERO'] = 0
    #
    #             # Write the reduced frame to disk
    #
    #             logger.debug(f"Saving processed image to {output_path}")
    #             proccessed_list.append(output_path)
    #             proc_hdu.writeto(output_path, overwrite=True)
    #
    #     return proccessed_list

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

    def export_observing_log(
            self,
            sub_dir: str = ""
    ):
        """Function to export observing log to file

        Parameters
        ----------
        sub_dir: subdirectory associated with data, e.g date

        Returns
        -------
        """
        log = self.load_observing_log(sub_dir=sub_dir)
        path = self.get_observing_log_path(sub_dir=sub_dir)
        logger.debug(f"Saving to {path}")
        log.to_csv(path)

    @staticmethod
    def get_observing_log_path(
            sub_dir: str = ""
    ):
        """Function to find path for observing log file

        Parameters
        ----------
        sub_dir: subdirectory associated with data, e.g date

        Returns
        -------
        path to the observing log
        """
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

        # Some fields should always go to the log

        key_list = self.header_keys + ["OBSCLASS"]

        for img_file in img_list:
            _, header = self.open_fits(img_file)

            row = []

            for key in key_list:
                row.append(header[key])

            row.append(img_file)

            log_data.append(row)

        log = pd.DataFrame(log_data, columns=key_list + [raw_img_key])
        log = log.sort_values(by=self.header_keys[0]).reset_index(drop=True)

        return log



