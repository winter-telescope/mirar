import logging
from abc import ABC

import astropy.io.fits
import numpy as np
import os
import socket
import getpass
import datetime
from winterdrp.io import save_to_path, open_fits
from winterdrp.paths import cal_output_sub_dir, get_mask_path, latest_save_key, latest_mask_save_key, get_output_path,\
    ProcessingError, base_name_key
from winterdrp.errors import ErrorReport

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

        self.night = None
        self.night_sub_dir = None
        self.preceding_steps = None

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
            batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]]
    ) -> tuple[list[list[list[np.ndarray], list[astropy.io.fits.header]]], list]:

        passed_batches = []
        failures = []

        for i, [images, headers] in enumerate(batches):

            try:
                images, headers = self._apply_to_images(images, headers)
                headers = self._update_processing_history(headers)
                passed_batches.append([images, headers])
            except ProcessingError as e:
                err = ErrorReport(e, self.__module__, images, headers)
                logger.error(err.generate_log_message())
                failures.append(err)

        batches = self.update_batches(passed_batches)

        return batches, failures

    @staticmethod
    def update_batches(
        batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]]
    ) -> list[list[list[np.ndarray], list[astropy.io.fits.header]]]:
        return batches

    def check_prerequisites(
            self,
    ):
        pass

    @staticmethod
    def open_fits(
            path: str
    ) -> tuple[np.ndarray, astropy.io.fits]:
        return open_fits(path)

    @staticmethod
    def save_fits(
            data,
            header,
            path,
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


class ProcessorWithCache(BaseProcessor, ABC):

    def __init__(
            self,
            try_load_cache: bool = False,
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

        return output_path

    def get_cache_file_name(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> str:
        raise NotImplementedError

    def get_cache_file(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[np.ndarray, astropy.io.fits.Header]:

        path = self.get_cache_path(images, headers)

        exists = os.path.exists(path)

        if np.logical_and(self.try_load_cache, exists):
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

# class ProcessorWithCache(BaseProcessor, ABC):
#
#     @staticmethod
#     def select_cache_images(
#             observing_log: pd.DataFrame
#     ) -> list:
#         raise NotImplementedError
#
#     def load_cache_file(
#             self,
#             path: str
#     ) -> (np.ndarray, astropy.io.fits.Header):
#
#         if path in self.cache:
#             img, header = self.cache[path]
#         else:
#             img, header = self.open_fits(path)
#             self.cache[path] = (img, header)
#         return img, header
#
#     def make_cache(
#             self,
#             observing_log: pd.DataFrame,
#     ):
#
#         img_path_list = self.select_cache_images(observing_log)
#
#         if len(img_path_list) > 0:
#             self.make_cache_files(
#                 img_path_list,
#             )
#
#     def make_cache_files(
#             self,
#             image_paths: list[str],
#     ):
#         raise NotImplementedError
#
#     def subselect_log(
#             self,
#             observing_log: pd.DataFrame,
#             key: str
#     ) -> pd.DataFrame:
#         mask = np.logical_and(
#             observing_log["TARGET"] == key,
#             observing_log["NIGHT"] == self.night
#         )
#
#         return observing_log[mask]
#
#     def select_from_log(
#             self,
#             observing_log: pd.DataFrame,
#             key: str
#     ) -> [str]:
#
#         obs = self.subselect_log(
#             observing_log=observing_log,
#             key=key
#         )
#
#         logger.debug(f"Found {len(obs)} entries with key '{key}' for night '{self.night}'")
#
#         return list(obs[raw_img_key])
