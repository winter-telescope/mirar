import logging
import os

import astropy.io.fits
import numpy as np
from astropy.io import fits
from glob import glob
import pandas as pd
import copy
from astropy.time import Time
from astropy import units as u
from winterdrp.paths import \
    cal_output_dir,\
    raw_img_dir, \
    observing_log_dir

logger = logging.getLogger(__name__)


raw_img_key = "RAWIMAGEPATH"
core_fields = ["OBSCLASS", "TARGET", "UTCTIME"]


class Pipeline:

    pipelines = {}
    name = None

    # Set up elements to use
    astrometry = ("GAIA", 9., 13.)
    photometry_cal = dict()

    @property
    def pipeline_configurations(self):
        raise NotImplementedError()

    # Fix keys from header to save in log
    # The first key should sort images in order

    header_keys = list()

    # Key from log which defines batch grouping
    # By default this is the raw image path, i.e each image is processed separately

    batch_split_keys = ["RAWIMAGEPATH"]

    standard_flats_dir = None

    default_log_history_nights = 0

    def __init__(
            self,
            pipeline_configuration: str = None,
            night: int | str = "",
            log_history_nights: int = None
    ):

        if log_history_nights is None:
            log_history_nights = self.default_log_history_nights

        self.night_sub_dir = os.path.join(self.name, night)

        self.processors = self.load_pipeline_configuration(pipeline_configuration)

        self.observing_logs_cache = dict()

        observing_logs = self.load_observing_log_block(
            night_sub_dir=self.night_sub_dir,
            log_history_nights=log_history_nights
        )

        self.configure_processors(sub_dir=self.night_sub_dir)

        for i, (processor) in enumerate(self.processors):
            # preceding_steps = [x[0] for x in self.processors[:i]]
            processor.set_preceding_steps(previous_steps=self.processors[:i])
            processor.make_cache(observing_log=observing_logs)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name in cls.pipelines.keys():
            err = f"Pipeline name '{cls.name}' is already found in the pipeline registered keys. " \
                  f"The Pipeline class variable 'name' must be unique!"
            logger.error(err)
            raise ValueError(err)
        cls.pipelines[cls.name] = cls

    def configure_processors(
            self,
            sub_dir: str = ""
    ):
        for processor in self.processors:
            processor.set_sub_dir(sub_dir=sub_dir)
            processor.set_open_fits(self.open_fits)

    def open_fits(
            self,
            path: str
    ) -> (np.array, astropy.io.fits.Header):
        with fits.open(path) as raw:
            data, header = self.reformat_raw_data(raw, path)

        for key in core_fields:
            if key not in header.keys():
                err = f"Essential key {key} not found in header. " \
                      f"Please add this field first. Available fields are: {header.keys()}"
                logger.error(err)
                raise KeyError(err)

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
            img: astropy.io.fits.HDUList,
            path: str
    ) -> (np.ndarray, astropy.io.fits.Header):
        raise NotImplementedError

    def load_pipeline_configuration(
            self,
            configuration: str | list = None,
    ):
        if isinstance(configuration, str | None):
            return copy.copy(self.pipeline_configurations[configuration])
        else:
            return copy.copy(configuration)

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

        for processor in self.processors:
            processor.make_cache(
                observing_log,
                sub_dir=sub_dir,
                preceding_steps=preceding_steps,
            )
            preceding_steps.append(processor.apply)

    def split_raw_images_into_batches(
            self,
            select_batch: str = None
    ) -> list:

        observing_log = self.load_observing_log(sub_dir=self.sub_dir)

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

        if select_batch is not None:
            if select_batch in batches:
                logger.debug(f"Only selecting '{select_batch}'")
                batches = [select_batch]
            else:
                err = f"The batch '{select_batch}' was selected, " \
                      f"but was not found in the batch list."
                logger.error(err)
                raise KeyError(err)

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

        for processor in self.processors:
            logger.debug(f"Applying '{processor}' processor to {len(images)} images")
            images, headers = processor.apply(
                images,
                headers,
                sub_dir=sub_dir
            )

        return images

    # def process_images(self, sub_dir="", raw_image_list=None, reprocess=True):
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

    def load_observing_log(
            self,
            sub_dir: str | int
    ):

        path = self.get_observing_log_path(sub_dir=sub_dir)

        if path in self.observing_logs_cache:
            log = self.observing_logs_cache[path]
        else:
            log = self.parse_observing_log(night_sub_dir=sub_dir)
            self.observing_logs_cache[path] = log
            self.export_observing_log(sub_dir=sub_dir)

        return log

    def load_observing_log_block(
            self,
            night_sub_dir: str | int,
            log_history_nights: int
    ) -> pd.DataFrame:

        log = self.load_observing_log(sub_dir=night_sub_dir)

        pipeline, night = night_sub_dir.split("/")

        if len(str(night)) != 8:
            err = f"Night format not recignised. Folders should be organised as YYYYMMDD. " \
                  f"Instead found {night}."
            logger.error(err)
            raise ValueError(err)

        date = Time(f"{night[:4]}-{night[4:6]}-{night[6:]}")

        other_dates = [
            (date - ((x+1) * u.day)).isot.split("T")[0].replace("-", "")
            for x in range(log_history_nights)[::-1]
        ]

        for other_date in other_dates:

            other_sub_dir = os.path.join(pipeline, other_date)

            try:
                log = pd.concat(
                    log,
                    self.load_observing_log(sub_dir=other_sub_dir)
                )
            except NotADirectoryError:
                msg = f"Did not find sub diretory {other_sub_dir}, skipping instead."
                logger.warning(msg)

        return log

    def parse_observing_log(self, night_sub_dir):

        raw_dir = raw_img_dir(night_sub_dir)

        if not os.path.isdir(raw_dir):
            error = f"Raw image directory '{raw_dir}' does not exit"
            logger.error(error)
            raise NotADirectoryError(error)

        img_list = glob(f'{raw_dir}/*.fits')

        if len(img_list) == 0:
            err = f"No images found in directory {raw_dir}"
            logger.error(err)
            raise FileNotFoundError(err)
        else:
            logger.debug(f"Found {len(img_list)} raw images in {raw_dir}")

        log_data = []

        # Some fields should always go to the log

        key_list = self.header_keys + core_fields + ["NIGHT"]

        for img_file in img_list:
            _, header = self.open_fits(img_file)

            header["NIGHT"] = night_sub_dir.split("/")[1]

            row = []

            for key in key_list:
                row.append(header[key])

            row.append(img_file)

            log_data.append(row)

        log = pd.DataFrame(log_data, columns=key_list + [raw_img_key])
        log = log.sort_values(by=self.header_keys[0]).reset_index(drop=True)

        return log
