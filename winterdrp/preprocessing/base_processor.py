import logging

import astropy.io.fits
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class BaseProcessor:

    base_name = None
    base_key = None

    def __init__(
            self,
            instrument_vars: dict,
            *args,
            **kwargs
    ):
        self.cache = dict()
        self.open_fits = instrument_vars["open_fits"]

    def apply_to_images(
            self,
            images: list,
            headers: list,
    ) -> list:
        raise NotImplementedError

    @staticmethod
    def get_file_path(header, sub_dir=""):
        raise NotImplementedError

    def load_cache_file(
            self,
            path: str
    ) -> astropy.io.fits.HDUList:

        if path in self.cache:
            img = self.cache[path]
        else:
            img = self.open_fits(path)
            self.cache[path] = img
        return img

    def select_cache_images(
            self,
            observing_log: pd.DataFrame
    ) -> list:
        mask = observing_log["OBJECT"] == self.base_key
        return observing_log[mask]["RAWIMAGEPATH"]

    def make_cache_files(
            self,
            image_list: list,
            preceding_steps: list,
            sub_dir: str = "",
            *args,
            **kwargs
    ):
        pass

    def save_fits(self, img, path):
        self.cache[path] = img.data
        logger.info(f"Saving to {path}")
        img.writeto(path, overwrite=True)
