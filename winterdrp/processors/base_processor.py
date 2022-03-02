import logging
from abc import ABC

import astropy.io.fits
import numpy as np
import pandas as pd
from winterdrp.io import create_fits

logger = logging.getLogger(__name__)


class BaseProcessor:

    base_name = None
    base_key = None

    requires = []

    subclasses = {}

    def __init__(
            self,
            instrument_vars: dict,
            *args,
            **kwargs
    ):
        self.cache = dict()
        self.open_fits = instrument_vars["open_fits"]

        # Check processor prerequisites are satisfied

        self.check_prerequisites(instrument_vars["preceding_steps"])

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def _apply_to_images(
            self,
            images: list,
            headers: list,
            sub_dir: str = ""
    ) -> (list, list):
        raise NotImplementedError

    def _update_processing_history(
            self,
            headers: list,
    ) -> list:
        for header in headers:
            header["CALSTEPS"] += self.base_key
        return headers

    def apply(
            self,
            images: list,
            headers: list,
            sub_dir: str = ""
    ) -> (list, list):
        images, headers = self._apply_to_images(images, headers, sub_dir=sub_dir)
        headers = self._update_processing_history(headers)
        return images, headers

    def check_prerequisites(
            self,
            preceding_steps: list,
    ):
        pass

    # @staticmethod
    # def get_file_path(header, sub_dir=""):
    #     raise NotImplementedError

    def make_cache(
            self,
            observing_log: pd.DataFrame,
            preceding_steps: list,
            sub_dir: str = "",
    ):
        pass

    def save_fits(self, data, header, path):
        self.cache[path] = (data, header)
        logger.info(f"Saving to {path}")
        img = create_fits(data, header=header)
        img.writeto(path, overwrite=True)


class ProcessorWithCache(BaseProcessor, ABC):

    @staticmethod
    def select_cache_images(
            observing_log: pd.DataFrame
    ) -> list:
        return []

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
            preceding_steps: list,
            sub_dir: str = "",
    ):
        image_list = self.select_cache_images(observing_log)

        if len(image_list) > 0:
            self.make_cache_files(
                image_list,
                preceding_steps=preceding_steps,
                sub_dir=sub_dir
            )

    def make_cache_files(
            self,
            image_list: list,
            preceding_steps: list,
            sub_dir: str = ""
    ):
        raise NotImplementedError
