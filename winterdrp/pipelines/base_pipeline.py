import logging
import os

import astropy.io.fits
import numpy as np
from glob import glob
import pandas as pd
import copy
from astropy.time import Time
from astropy import units as u
from winterdrp.paths import cal_output_dir, raw_img_dir, observing_log_dir, raw_img_key, saturate_key, \
    get_preprocess_path, ProcessingError
from winterdrp.processors.base_processor import ProcessorWithCache
from winterdrp.io import save_to_path


logger = logging.getLogger(__name__)


core_fields = ["OBSCLASS", "TARGET", "UTCTIME"]


class Pipeline:

    pipelines = {}
    name = None

    @property
    def pipeline_configurations(self):
        raise NotImplementedError()

    @property
    def gain(self):
        raise NotImplementedError()

    @property
    def non_linear_level(self):
        raise NotImplementedError()

    # Fix keys from header to save in log
    # The first key should sort images in order

    header_keys = list()

    # Key from log which defines batch grouping
    # By default this is the raw image path, i.e each image is processed separately

    batch_split_keys = ["RAWIMAGEPATH"]

    default_log_history_nights = 0

    def __init__(
            self,
            pipeline_configuration: str = None,
            night: int | str = "",
            log_history_nights: int = None,
            skip_build_cache: bool = False,
            remake_logs: bool = False
    ):
        self.remake_logs = remake_logs

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

        if skip_build_cache:
            logger.warning("Skipping cache building for all processors.")

        for i, (processor) in enumerate(self.processors):

            logger.debug(f"Initialising processor {processor.__class__}")
            processor.set_preceding_steps(previous_steps=self.processors[:i])
            processor.check_prerequisites()

        logger.debug("Pipeline initialisation complete.")

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
            processor.set_night(night_sub_dir=sub_dir)

    def open_raw_image(
            self,
            path: str
    ) -> tuple[np.array, astropy.io.fits.Header]:

        data, header = self.load_raw_image(path)

        for key in core_fields:
            if key not in header.keys():
                err = f"Essential key {key} not found in header. " \
                      f"Please add this field first. Available fields are: {list(header.keys())}"
                logger.error(err)
                raise KeyError(err)

        return data.astype(np.float64), header

    def open_raw_image_batch(
            self,
            paths: list
    ) -> tuple[list, list]:

        images = []
        headers = []
        for path in paths:
            data, header = self.open_raw_image(path)
            images.append(data)
            headers.append(header)

        return images, headers

    @staticmethod
    def load_raw_image(
            path: str
    ) -> tuple[np.array, astropy.io.fits.Header]:
        raise NotImplementedError

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        raise NotImplemented

    def load_pipeline_configuration(
            self,
            configuration: str | list = None,
    ):
        if isinstance(configuration, str | None):
            return copy.copy(self.pipeline_configurations[configuration])
        else:
            return copy.copy(configuration)

    # def make_calibration_files(self, sub_dir=""):
    #
    #     observing_log = self.load_observing_log(night_sub_dir=sub_dir)
    #
    #     cal_dir = cal_output_dir(sub_dir)
    #
    #     # Make calibration directory, unless it already exists
    #
    #     try:
    #         os.makedirs(cal_dir)
    #     except OSError:
    #         pass
    #
    #     logger.info(f"Making calibration files for directory {raw_img_dir(sub_dir)}")
    #
    #     preceding_steps = []
    #
    #     for processor in self.processors:
    #         processor.make_cache(
    #             observing_log,
    #             sub_dir=sub_dir,
    #             preceding_steps=preceding_steps,
    #         )
    #         preceding_steps.append(processor.apply)

    def split_raw_images_into_batches(
            self,
            select_batch: str = None
    ) -> list:

        observing_log = self.load_observing_log(night_sub_dir=self.night_sub_dir)

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
            batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]],
    ):

        for processor in self.processors:
            logger.debug(f"Applying '{processor.__class__}' processor to {len(batches)} batches")
            batches, failures = processor.apply(
                batches
            )

        return batches

    def export_observing_log(
            self,
            night_sub_dir: str = ""
    ):
        """Function to export observing log to file

        Parameters
        ----------
        night_sub_dir: subdirectory associated with data, e.g date

        Returns
        -------
        """
        log = self.load_observing_log(night_sub_dir=night_sub_dir)
        path = self.get_observing_log_path(night_sub_dir=night_sub_dir)
        logger.debug(f"Saving to {path}")
        log.to_csv(path)

    @staticmethod
    def get_observing_log_path(
            night_sub_dir: str | int = ""
    ):
        """Function to find path for observing log file

        Parameters
        ----------
        night_sub_dir: subdirectory associated with data, e.g date

        Returns
        -------
        path to the observing log
        """
        base_dir = observing_log_dir(sub_dir=night_sub_dir)
        return os.path.join(base_dir, "observing_log.csv")

    def load_observing_log(
            self,
            night_sub_dir: str | int
    ):

        path = self.get_observing_log_path(night_sub_dir=night_sub_dir)

        if path in self.observing_logs_cache:
            log = self.observing_logs_cache[path]
        elif np.logical_and(not self.remake_logs, os.path.isfile(path)):
            logger.debug(f"Found log, loading {path}")
            log = pd.read_csv(path)
            self.observing_logs_cache[path] = log
        else:
            logger.debug(f"Making log {path}")
            log = self.parse_observing_log(night_sub_dir=night_sub_dir)
            self.observing_logs_cache[path] = log
            self.export_observing_log(night_sub_dir=night_sub_dir)

        return log

    def load_observing_log_block(
            self,
            night_sub_dir: str | int,
            log_history_nights: int
    ) -> pd.DataFrame:

        log = self.load_observing_log(night_sub_dir=night_sub_dir)

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
                    self.load_observing_log(night_sub_dir=other_sub_dir)
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

        preprocess_dir = os.path.dirname(get_preprocess_path(img_list[0]))

        try:
            os.makedirs(preprocess_dir)
        except OSError:
            pass

        log_data = []

        # Some fields should always go to the log

        key_list = self.header_keys + core_fields + ["NIGHT"]

        for raw_img_path in img_list:
            data, header = self.open_raw_image(raw_img_path)

            preprocess_img_path = get_preprocess_path(raw_img_path)

            header["NIGHT"] = night_sub_dir.split("/")[1]

            save_to_path(data, header, preprocess_img_path)

            row = []

            for key in key_list:
                row.append(header[key])

            row.append(preprocess_img_path)

            log_data.append(row)

        log = pd.DataFrame(log_data, columns=key_list + [raw_img_key])
        log = log.sort_values(by=self.header_keys[0]).reset_index(drop=True)

        return log

    def set_saturation(
            self,
            header: astropy.io.fits.Header
    ) -> astropy.io.fits.Header:
        # update the SATURATE keyword in the header for subsequent sextractor runs
        co_add_head = header['COADDS']
        num_co_adds = int(co_add_head)
        saturation_level = self.non_linear_level * num_co_adds
        if "SKMEDSUB" in header.keys():
            saturation_level -= header['SKMEDSUB']
        header.append((saturate_key, saturation_level, 'Saturation level'), end=True)
        return header
