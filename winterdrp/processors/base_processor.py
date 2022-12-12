"""
Module containing the :class:`~wintedrp.processors.BaseProcessor`
"""
import datetime
import getpass
import hashlib
import logging
import os
import socket
from abc import ABC
from pathlib import Path
from queue import Queue
from threading import Thread

import numpy as np

from winterdrp.data import DataBatch, Dataset, Image, ImageBatch, SourceBatch
from winterdrp.errors import (
    ErrorReport,
    ErrorStack,
    NoncriticalProcessingError,
    ProcessorError,
)
from winterdrp.io import open_fits, save_to_path
from winterdrp.paths import (
    base_name_key,
    cal_output_sub_dir,
    get_mask_path,
    get_output_path,
    latest_mask_save_key,
    latest_save_key,
    max_n_cpu,
    package_name,
    proc_history_key,
    raw_img_key,
)

logger = logging.getLogger(__name__)


class PrerequisiteError(ProcessorError):
    """
    An error raised if a processor requires another one as a prerequisite,
    but that processor is not present
    """


class NoCandidatesError(ProcessorError):
    """
    An error raised if a :class:`~wintedrp.processors.CandidateGenerator` produces
    no candidates
    """


class BaseProcessor:
    """
    Base processor class, to be inherited from for all processors
    """

    @property
    def base_key(self):
        """
        Unique key for the processor, to be used e.g in processing history tracking

        :return: None
        """
        raise NotImplementedError

    max_n_cpu: int = max_n_cpu

    subclasses = {}

    def __init__(self):
        self.night = None
        self.night_sub_dir = None
        self.preceding_steps = None
        self.passed_dataset = self.err_stack = None

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def set_preceding_steps(self, previous_steps: list):
        """
        Provides processor with the list of preceding processors, and saves this

        :param previous_steps: list of processors
        :return: None
        """
        self.preceding_steps = previous_steps

    def set_night(self, night_sub_dir: str | int = ""):
        """
        Sets the night subdirectory for the processor to read/write data

        :param night_sub_dir: String/int subdirectory for night
        :return: None
        """
        self.night_sub_dir = night_sub_dir
        self.night = night_sub_dir.split("/")[-1]

    def generate_error_report(
        self, exception: Exception, batch: DataBatch
    ) -> ErrorReport:
        return ErrorReport(exception, self.__module__, batch.get_raw_image_names())

    def update_dataset(self, dataset: Dataset) -> Dataset:
        return dataset

    def check_prerequisites(
        self,
    ):
        pass

    def clean_cache(self):
        """
        Function to clean the internal cache filled by base_apply

        :return: None
        """
        self.passed_dataset = self.err_stack = None

    def base_apply(self, dataset: Dataset) -> tuple[Dataset, ErrorStack]:
        """
        Core function to act on a dataset, and return an updated dataset

        :param dataset: Input dataset
        :return: Updated dataset, and any caught errors
        """

        self.passed_dataset = Dataset()
        self.err_stack = ErrorStack()

        if len(dataset) > 0:

            n_cpu = min([self.max_n_cpu, len(dataset)])

            watchdog_queue = Queue()

            workers = []

            for _ in range(n_cpu):
                # Set up a worker thread to process database load
                worker = Thread(target=self.apply_to_batch, args=(watchdog_queue,))
                worker.daemon = True
                worker.start()

                workers.append(worker)

            for batch in dataset:
                watchdog_queue.put(item=batch)

            watchdog_queue.join()

        dataset = self.update_dataset(self.passed_dataset)
        err_stack = self.err_stack

        self.clean_cache()

        return dataset, err_stack

    def apply_to_batch(self, queue):
        """
        Function to run self.apply on a batch in the queue, catch any errors, and then
        update the internal cache with the results.

        :param queue: python threading queue
        :return: None
        """
        while True:
            batch = queue.get()
            try:
                batch = self.apply(batch)
                self.passed_dataset.append(batch)
            except NoncriticalProcessingError as exc:
                err = self.generate_error_report(exc, batch)
                logger.error(err.generate_log_message())
                self.err_stack.add_report(err)
                self.passed_dataset.append(batch)
            except Exception as exc:  # pylint: disable=broad-except
                err = self.generate_error_report(exc, batch)
                logger.error(err.generate_log_message())
                self.err_stack.add_report(err)
            queue.task_done()

    def apply(self, batch: DataBatch):
        """
        Function applying the processor to a
        :class:`~winterdrp.data.base_data.DataBatch`.
        Also updates the processing history.

        :param batch: input data batch
        :return: updated data batch
        """
        batch = self._apply(batch)
        batch = self._update_processing_history(batch)
        return batch

    def _apply(self, batch: DataBatch) -> DataBatch:
        """
        Core function to update the :class:`~winterdrp.data.base_data.DataBatch`

        :param batch: Input data batch
        :return: updated data batch
        """
        raise NotImplementedError

    def _update_processing_history(
        self,
        batch: DataBatch,
    ) -> DataBatch:
        """
        Function to update the processing history of each
        :class:`~winterdrp.data.base_data.DataBlock` object in a
        :class:`~winterdrp.data.base_data.DataBatch`.

        :param batch: Input data batch
        :return: Updated data batch
        """
        for i, data_block in enumerate(batch):
            data_block[proc_history_key] += self.base_key + ","
            data_block["REDUCER"] = getpass.getuser()
            data_block["REDMACH"] = socket.gethostname()
            data_block["REDTIME"] = str(datetime.datetime.now())
            data_block["REDSOFT"] = package_name
            batch[i] = data_block
        return batch


class CleanupProcessor(BaseProcessor, ABC):
    def update_dataset(self, dataset: Dataset) -> Dataset:
        # Remove empty dataset
        new_dataset = Dataset([x for x in dataset.get_batches() if len(x) > 0])
        return new_dataset


class ImageHandler:
    @staticmethod
    def open_fits(path: str | Path) -> Image:
        path = str(path)
        data, header = open_fits(path)
        if raw_img_key not in header:
            header[raw_img_key] = path
        if base_name_key not in header:
            header[base_name_key] = Path(path).name
        return Image(data=data, header=header)

    @staticmethod
    def save_fits(
        image: Image,
        path: str | Path,
    ):
        path = str(path)
        data = image.get_data()
        header = image.get_header()
        if header is not None:
            header[latest_save_key] = path
        logger.info(f"Saving to {path}")
        save_to_path(data, header, path)

    def save_mask(self, image: Image, img_path: str) -> str:
        data = image.get_data()
        mask = (~np.isnan(data)).astype(float)
        mask_path = get_mask_path(img_path)
        header = image.get_header()
        header[latest_mask_save_key] = mask_path
        self.save_fits(Image(mask, header), mask_path)
        return mask_path

    @staticmethod
    def get_hash(image: ImageBatch):
        key = "".join(sorted([x[base_name_key] + x[proc_history_key] for x in image]))
        return hashlib.sha1(key.encode()).hexdigest()


class BaseImageProcessor(BaseProcessor, ImageHandler, ABC):
    def _apply(self, batch: ImageBatch) -> ImageBatch:
        return self._apply_to_images(batch)

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        raise NotImplementedError


class ProcessorWithCache(BaseImageProcessor, ABC):
    def __init__(
        self,
        *args,
        try_load_cache: bool = True,
        write_to_cache: bool = True,
        overwrite: bool = True,
        cache_sub_dir: str = cal_output_sub_dir,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.try_load_cache = try_load_cache
        self.write_to_cache = write_to_cache
        self.overwrite = overwrite
        self.cache_sub_dir = cache_sub_dir

    def select_cache_images(self, images: ImageBatch) -> ImageBatch:
        raise NotImplementedError

    def get_cache_path(self, images: ImageBatch) -> str:

        file_name = self.get_cache_file_name(images)

        output_path = get_output_path(
            base_name=file_name, dir_root=self.cache_sub_dir, sub_dir=self.night_sub_dir
        )

        try:
            os.makedirs(os.path.dirname(output_path))
        except OSError:
            pass

        return output_path

    def get_cache_file_name(self, images: ImageBatch) -> str:
        cache_images = self.select_cache_images(images)
        return f"{self.base_key}_{self.get_hash(cache_images)}.fits"

    def get_cache_file(self, images: ImageBatch) -> Image:

        path = self.get_cache_path(images)

        exists = os.path.exists(path)

        if np.logical_and(self.try_load_cache, exists):
            logger.info(f"Loading cached file {path}")
            return self.open_fits(path)

        image = self.make_image(images)

        if self.write_to_cache:
            if np.sum([not exists, self.overwrite]) > 0:
                self.save_fits(image, path)

        return image

    def make_image(self, images: ImageBatch) -> Image:
        raise NotImplementedError


class ProcessorPremadeCache(ProcessorWithCache, ABC):
    def __init__(self, master_image_path: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.master_image_path = master_image_path

    def get_cache_path(self, images: ImageBatch) -> str:
        return self.master_image_path


class BaseCandidateGenerator(BaseProcessor, ImageHandler, ABC):
    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def _apply(self, batch: ImageBatch) -> SourceBatch:
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

    def _apply(self, batch: SourceBatch) -> SourceBatch:
        return self._apply_to_candidates(batch)

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        raise NotImplementedError
