import logging
from abc import ABC
import astropy.io.fits
import numpy as np
import os
import socket
import getpass
import datetime
import hashlib
from pathlib import Path

from winterdrp.io import save_to_path, open_fits
from winterdrp.paths import cal_output_sub_dir, get_mask_path, latest_save_key, latest_mask_save_key, get_output_path,\
    base_name_key, proc_history_key, raw_img_key
from winterdrp.errors import ErrorReport, ErrorStack, ProcessorError, NoncriticalProcessingError
from winterdrp.data import DataBatch, DataSet, Image, ImageBatch, SourceBatch

logger = logging.getLogger(__name__)


class PrerequisiteError(ProcessorError):
    pass


class NoCandidatesError(ProcessorError):
    pass


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
        batches: list[DataBatch]
    ) -> list[DataBatch]:
        return batches

    def check_prerequisites(
            self,
    ):
        pass

    def base_apply(
            self,
            dataset: DataSet
    ) -> tuple[dataset, ErrorStack]:

        passed_batches = []
        err_stack = ErrorStack()

        for i, batch in enumerate(dataset):

            try:
                batch = self.apply(batch)
                passed_batches.append(batch)
            except NoncriticalProcessingError as e:
                err = self.generate_error_report(e, batch)
                logger.error(err.generate_log_message())
                err_stack.add_report(err)
                passed_batches.append(batch)
            except Exception as e:
                err = self.generate_error_report(e, batch)
                logger.error(err.generate_log_message())
                err_stack.add_report(err)

        batches = self.update_batches(passed_batches)

        return batches, err_stack

    def apply(self, batch: DataBatch):
        raise NotImplementedError

    def generate_error_report(self, exception: Exception, batch: DataBatch) -> ErrorReport:
        raise NotImplementedError


class ImageHandler:

    @staticmethod
    def open_fits(
            path: str
    ) -> Image:
        image, header = open_fits(path)
        return Image(image=image, header=header)

    @staticmethod
    def save_fits(
            image,
            path: str,
    ):
        data = image.get_data()
        header = image.get_header()
        if header is not None:
            header[latest_save_key] = path
        logger.info(f"Saving to {path}")
        save_to_path(data, header, path)

    def save_mask(
            self,
            image: Image,
            img_path: str
    ) -> str:
        data = image.get_data()
        mask = (~np.isnan(data)).astype(float)
        mask_path = get_mask_path(img_path)
        header = image.get_header()
        header[latest_mask_save_key] = mask_path
        self.save_fits(Image(mask, header), mask_path)
        return mask_path

    @staticmethod
    def get_hash(headers: list[astropy.io.fits.Header]):
        key = "".join(sorted([x[base_name_key] + x[proc_history_key] for x in headers]))
        return hashlib.sha1(key.encode()).hexdigest()

    def image_batch_error_report(self, exception: Exception, batch: ImageBatch):
        contents = [Path(x).name for x in ",".join([x[raw_img_key] for x in batch[1]]).split(",")]
        return ErrorReport(exception, self.__module__, contents)


class BaseImageProcessor(BaseProcessor, ImageHandler, ABC):

    def apply(
            self,
            batch: ImageBatch
    ) -> ImageBatch:
        batch = self._apply_to_images(batch)
        batch = self._update_processing_history(batch)

        return batch

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:
        raise NotImplementedError

    def _update_processing_history(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:
        for i, image in enumerate(batch):
            image[proc_history_key] += self.base_key + ","
            image['REDUCER'] = getpass.getuser()
            image['REDMACH'] = socket.gethostname()
            image['REDTIME'] = str(datetime.datetime.now())
            image["REDSOFT"] = "winterdrp"
            batch[i] = image
        return batch

    def generate_error_report(self, exception: Exception, batch: ImageBatch) -> ErrorReport:
        return self.image_batch_error_report(exception, batch)


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

    def select_cache_images(self, images: list[Image]) -> list[Image]:
        raise NotImplementedError

    def get_cache_path(self, images: list[Image]) -> str:

        file_name = self.get_cache_file_name(images)

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

    def get_cache_file_name(self, images: list[Image]) -> str:
        cache_images = self.select_cache_images(images)
        return f"{self.base_key}_{self.get_hash(cache_images)}.fits"

    def get_cache_file(self, images: list[Image]) -> Image:

        path = self.get_cache_path(images)

        exists = os.path.exists(path)

        if np.logical_and(self.try_load_cache, exists):
            logger.info(f"Loading cached file {path}")
            return self.open_fits(path)

        else:

            image = self.make_image(images)

            if self.write_to_cache:
                if np.sum([not exists, self.overwrite]) > 0:
                    self.save_fits(image, path)

        return image

    def make_image(self, images: list[Image]) -> Image:
        raise NotImplementedError


class ProcessorPremadeCache(ProcessorWithCache, ABC):

    def __init__(
            self,
            master_image_path: str,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.master_image_path = master_image_path

    def get_cache_path(self, images: list[Image]) -> str:
        return self.master_image_path


class BaseCandidateGenerator(BaseProcessor, ImageHandler, ABC):

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def apply(self, batch: ImageBatch) -> SourceBatch:
        source_batch = self._apply_to_images(batch)

        if len(source_batch) == 0:
            err = "No sources found in image batch"
            logger.error(err)
            raise NoCandidatesError(err)

        return source_batch

    def _apply_to_images(self, batch: ImageBatch) -> SourceBatch:
        raise NotImplementedError


class BaseDataframeProcessor(BaseProcessor, ABC):

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def apply(
            self,
            batch: SourceBatch
    ) -> SourceBatch:
        return batch

    def _apply_to_candidates(
            self,
            source_list: SourceBatch,
    ) -> SourceBatch:
        raise NotImplementedError

    def generate_error_report(self, exception: Exception, batch: SourceBatch):
        contents = batch[base_name_key]
        return ErrorReport(exception, self.__module__, contents)
