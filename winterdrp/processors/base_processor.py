import logging
from abc import ABC
from typing import Tuple

import astropy.io.fits
import numpy as np
import os
import socket
import getpass
import datetime
import hashlib
import pandas as pd

from winterdrp.io import save_to_path, open_fits
from winterdrp.paths import cal_output_sub_dir, get_mask_path, latest_save_key, latest_mask_save_key, get_output_path,\
    ProcessingError, base_name_key, proc_history_key
from winterdrp.errors import ErrorReport

logger = logging.getLogger(__name__)


class PrerequisiteError(BaseException):
    pass

class ImageHandler:
    @staticmethod
    def open_fits(
            path: str
    ) -> tuple[np.ndarray, astropy.io.fits]:
        return open_fits(path)

    @staticmethod
    def save_fits(
            data,
            header,
            path: str,
    ):
        header[latest_save_key] = path
        logger.info(f"Saving to {path}")
        save_to_path(data, header, path)

    def save_mask(
            self,
            data: np.ndarray,
            header: astropy.io.fits.Header,
            img_path: str
    ) -> str:
        mask = (~np.isnan(data)).astype(float)
        mask_path = get_mask_path(img_path)
        header[latest_mask_save_key] = mask_path
        self.save_fits(mask, header, mask_path)
        return mask_path

    @staticmethod
    def get_hash(headers: list[astropy.io.fits.Header]):
        key = "".join(sorted([x[base_name_key] + x[proc_history_key] for x in headers]))
        return hashlib.sha1(key.encode()).hexdigest()


class BaseProcessor:

    @property
    def base_key(self):
        raise NotImplementedError()

    subclasses = {}

    def __init__(
            self,
            *args,
            **kwargs
    ):

        self.night = None
        self.night_sub_dir = None
        self.preceding_steps = None

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def set_preceding_steps(
            self,
            previous_steps: list
    ):
        self.preceding_steps = previous_steps

    def set_night(
            self,
            night_sub_dir: str | int = ""
    ):
        self.night_sub_dir = night_sub_dir
        self.night = night_sub_dir.split("/")[-1]

    @staticmethod
    def update_batches(
        batches: list
    ) -> list:
        return batches

    def check_prerequisites(
            self,
    ):
        pass

    def base_apply(
            self,
            batches: list
    ) -> tuple[list, list[ErrorReport]]:

        passed_batches = []
        failures = []

        for i, batch in enumerate(batches):

            try:
                batch = self.apply(batch)
                passed_batches.append(batch)
            except ProcessingError as e:
                err = ErrorReport(e, self.__module__, batch)
                logger.error(err.generate_log_message())
                failures.append(err)

        batches = self.update_batches(passed_batches)

        return batches, failures

    def apply(self, batch):
        raise NotImplementedError


class BaseImageProcessor(BaseProcessor, ImageHandler, ABC):

    def apply(
            self,
            batch: list[list[np.ndarray], list[astropy.io.fits.header]]
    ) -> list[list[np.ndarray], list[astropy.io.fits.header]]:

        [images, headers] = batch

        images, headers = self._apply_to_images(images, headers)
        headers = self._update_processing_history(headers)

        return [images, headers]

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        raise NotImplementedError

    def _update_processing_history(
            self,
            headers: list,
    ) -> list:
        for header in headers:
            header[proc_history_key] += self.base_key + ","
            header['REDUCER'] = getpass.getuser()
            header['REDMACH'] = socket.gethostname()
            header['REDTIME'] = str(datetime.datetime.now())
            header["REDSOFT"] = "winterdrp"
        return headers


class ProcessorWithCache(BaseImageProcessor, ABC):

    def __init__(
            self,
            try_load_cache: bool = True,
            write_to_cache: bool = True,
            overwrite: bool = True,
            cache_sub_dir: str = cal_output_sub_dir,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.try_load_cache = try_load_cache
        self.write_to_cache = write_to_cache
        self.overwrite = overwrite
        self.cache_sub_dir = cache_sub_dir

    def get_cache_path(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> str:

        file_name = self.get_cache_file_name(images, headers)

        output_path = get_output_path(
            base_name=file_name,
            dir_root=self.cache_sub_dir,
            sub_dir=self.night_sub_dir
        )

        try:
            os.makedirs(os.path.dirname(output_path))
        except OSError:
            pass

        return output_path

    def get_cache_file_name(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> str:
        return f"{self.base_key}_{self.get_hash(headers)}.fits"

    def get_cache_file(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[np.ndarray, astropy.io.fits.Header]:

        path = self.get_cache_path(images, headers)

        exists = os.path.exists(path)

        if np.logical_and(self.try_load_cache, exists):
            logger.info(f"Loading cached file {path}")
            return self.open_fits(path)

        else:

            image, header = self.make_image(images, headers)

            if self.write_to_cache:
                if np.sum([not exists, self.overwrite]) > 0:
                    self.save_fits(image, header, path)

        return image, header

    def make_image(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[np.ndarray, astropy.io.fits.Header]:
        raise NotImplementedError


class BaseCandidateGenerator(BaseProcessor, ImageHandler, ABC):

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def apply(
            self,
            batch: tuple[list[np.ndarray], list[astropy.io.fits.Header]]
    ) -> pd.DataFrame:
        [images, headers] = batch
        candidate_table = self._apply_to_images(images, headers)
        return candidate_table

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> pd.DataFrame:
        raise NotImplementedError


class BaseDataframeProcessor(BaseProcessor, ABC):

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def apply(
            self,
            batch: pd.DataFrame
    ) -> pd.DataFrame:
        batch = self._apply_to_candidates(batch)
        return batch

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        raise NotImplementedError
