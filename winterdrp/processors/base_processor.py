import logging
from abc import ABC

import astropy.io.fits
import numpy as np
import pandas as pd
import socket
import getpass
import datetime
from winterdrp.io import create_fits

logger = logging.getLogger(__name__)


class BaseProcessor:

    @property
    def base_key(self):
        raise NotImplementedError()

    requires = []

    subclasses = {}

    def __init__(
            self,
            *args,
            **kwargs
    ):
        self.cache = dict()

        self.night = None
        self.night_sub_dir = None
        self.open_fits = None
        self.preceding_steps = None

        # Check processor prerequisites are satisfied

        self.check_prerequisites()

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        raise NotImplementedError

    # def apply_previous(
    #         self,
    #         images: list[np.ndarray],
    #         headers: list[astropy.io.fits.header],
    # ) -> (list[np.ndarray], list[astropy.io.fits.header]):
    #     raise NotImplementedError

    def set_preceding_steps(
            self,
            previous_steps: list
    ):
        self.preceding_steps = previous_steps

    def set_open_fits(
            self,
            open_fits
    ):
        self.open_fits = open_fits

    def set_night(
            self,
            night_sub_dir: str | int = ""
    ):
        self.night_sub_dir = night_sub_dir
        self.night = night_sub_dir.split("/")[-1]

    def _update_processing_history(
            self,
            headers: list,
    ) -> list:
        for header in headers:
            header["CALSTEPS"] += self.base_key + ","
            header['REDUCER'] = getpass.getuser()
            header['REDMACH'] = socket.gethostname()
            header['REDTIME'] = str(datetime.datetime.now())
            header["REDSOFT"] = "winterdrp"
        return headers

    def apply(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        images, headers = self._apply_to_images(images, headers)
        headers = self._update_processing_history(headers)
        return images, headers

    def check_prerequisites(
            self,
    ):
        pass

    def make_cache(
            self,
            observing_log: pd.DataFrame,
    ):
        pass

    def save_fits(self, data, header, path):
        self.cache[path] = (data, header)
        logger.info(f"Saving to {path}")
        img = create_fits(data, header=header)
        img.writeto(path, overwrite=True)

    def load_and_apply_previous(
            self,
            img_path: str
    ) -> tuple[np.ndarray, astropy.io.fits.Header]:

        img, header = self.open_fits(img_path)

        # Iteratively apply corrections
        for p in self.preceding_steps:
            [img], [header] = p.apply([img], [header])

        return img, header


class ProcessorWithCache(BaseProcessor, ABC):

    @staticmethod
    def select_cache_images(
            observing_log: pd.DataFrame
    ) -> list:
        raise NotImplementedError

    def load_cache_file(
            self,
            path: str
    ) -> (np.ndarray, astropy.io.fits.Header):

        if path in self.cache:
            img, header = self.cache[path]
        else:
            img, header = self.open_fits(path)
            self.cache[path] = (img, header)
        return img, header

    def make_cache(
            self,
            observing_log: pd.DataFrame,
    ):

        img_path_list = self.select_cache_images(observing_log)

        if len(img_path_list) > 0:
            self.make_cache_files(
                img_path_list,
            )

    def make_cache_files(
            self,
            image_paths: list[str],
    ):
        raise NotImplementedError

    def subselect_log(
            self,
            observing_log: pd.DataFrame,
            key: str
    ) -> pd.DataFrame:
        mask = np.logical_and(
            observing_log["TARGET"] == key,
            observing_log["NIGHT"] == self.night
        )

        return observing_log[mask]

    def select_from_log(
            self,
            observing_log: pd.DataFrame,
            key: str
    ) -> [str]:

        obs = self.subselect_log(
            observing_log=observing_log,
            key=key
        )

        logger.debug(f"Found {len(obs)} entries with key '{key}' for night '{self.night}'")

        return list(obs["RAWIMAGEPATH"])
