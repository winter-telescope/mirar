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

    pipelines = {}
    name = None

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

        instrument_vars = dict([(x, getattr(self, x)) for x in dir(self)])

        self.observing_logs_cache = dict()
        self.processors = dict()

        for step in self.image_steps:

            if isinstance(step, str):
                name = step
            else:
                name = step[0]

                if len(step) > 1:
                    print(step[1:])

            print(step)
            #
            # if isinstance(processor_args, tuple):
            #     name = processor_args[0]
            #     args += processor_args[1:]
            # else:
            #     name = processor_args
            self.processors[name] = get_processor(name, instrument_vars, *args, **kwargs)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name in cls.pipelines.keys():
            err = f"Pipeline name '{cls.name}' is already found in the pipeline registered keys. " \
                  f"The Pipeline class variable 'name' must be unique!"
            logger.error(err)
            raise ValueError(err)
        cls.pipelines[cls.name] = cls

    def open_fits(
            self,
            path: str
    ) -> (np.array, astropy.io.fits.Header):
        with fits.open(path) as raw:
            data, header = self.reformat_raw_data(raw)
        return data, header

    def open_image_batch(
            self,
            paths: list
    ) -> (list, list):

        images = []
        headers = []
        for path in paths:
            data, header = self.open_fits(path)
            images.append(data)
            headers.append(header)

        return images, headers

    @staticmethod
    def reformat_raw_data(
            img: astropy.io.fits.HDUList
    ) -> (np.ndarray, astropy.io.fits.Header):
        raise NotImplementedError

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
            preceding_steps.append(self.processors[processor_name].apply)

    def split_raw_images_into_batches(
            self,
            sub_dir: str = ""
    ) -> list:
        observing_log = self.load_observing_log(sub_dir=sub_dir)

        mask = observing_log["OBSCLASS"] == "science"
        obs = observing_log[mask]

        split_vals = []

        for row in obs.itertuples():
            sv = []
            for key in self.batch_split_keys:
                sv.append(str(getattr(row, key)))

            split_vals.append("_".join(sv))

        split_vals = np.array(split_vals)

        batches = sorted(list(set(split_vals)))

        logger.debug(f"Selecting unique combinations of {self.batch_split_keys}, "
                     f"found the following batches: {batches}")

        raw_image_path_batches = [
            list(obs[split_vals == x][raw_img_key])
            for x in batches
        ]

        return sorted(raw_image_path_batches)

    def reduce_images(
            self,
            images: list,
            headers: list,
            sub_dir: str = ""
    ) -> list:

        for processor in self.image_steps:
            logger.debug(f"Applying '{processor}' processor to {len(images)} images")
            images, headers = self.processors[processor].apply(
                images,
                headers,
                sub_dir=sub_dir
            )

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

    # @staticmethod
    # def apply_astrometry(sub_dir="", redux_image_list=None, reprocess=True):
    #
    #     if redux_image_list is None:
    #         redux_image_list = parse_image_list(sub_dir, group_by_object=False, base_dir_f=reduced_img_dir)
    #
    #     # Try making output directory, unless it exists
    #
    #     output_dir = astrometry_output_dir(sub_dir)
    #
    #     try:
    #         os.makedirs(output_dir)
    #     except OSError:
    #         pass
    #
    #     # First run Sextractor
    #
    #     run_sextractor(
    #         redux_image_list,
    #         output_dir=output_dir,
    #         reprocess=reprocess
    #     )

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



