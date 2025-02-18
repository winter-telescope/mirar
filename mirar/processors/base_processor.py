"""
Module containing the :class:`~wintedrp.processors.BaseProcessor`
"""

import datetime
import getpass
import hashlib
import logging
import socket
import threading
from abc import ABC
from pathlib import Path
from queue import Queue
from threading import Thread
from typing import Callable

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from mirar.data import DataBatch, Dataset, Image, ImageBatch, SourceBatch
from mirar.errors import (
    ErrorReport,
    ErrorStack,
    NoncriticalProcessingError,
    ProcessorError,
)
from mirar.io import MissingCoreFieldError, open_fits, save_fits
from mirar.paths import (
    BASE_NAME_KEY,
    CAL_OUTPUT_SUB_DIR,
    LATEST_WEIGHT_SAVE_KEY,
    PACKAGE_NAME,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    core_source_fields,
    get_mask_path,
    get_output_path,
    max_n_cpu,
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


class DuplicationError(ProcessorError):
    """
    An error raised if a processor returns multiple images with the same name
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

        # For caching/multithreading
        self.passed_batches = {}
        self.err_stack = {}
        self.progress = {}

        # For tracking processing history
        self.latest_n_input_blocks = 0
        self.latest_n_input_batches = 0
        self.latest_n_output_blocks = 0
        self.latest_n_output_batches = 0
        self.processor_run_time_minutes = 0
        self.latest_error_stack = ErrorStack()

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def description(self) -> str:
        """
        Return a description of the processor

        :return: A description of the processor
        """
        raise NotImplementedError

    def __str__(self) -> str:
        return f"[{self.description()}]"

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
        """
        Generates an error report based on a python Exception

        :param exception: exception raised
        :param batch: batch which generated exception
        :return: error report
        """
        return ErrorReport(exception, self.__module__, batch.get_raw_image_names())

    def update_dataset(self, dataset: Dataset) -> Dataset:
        """
        Update a dataset after processing

        :param dataset: Initial dataset
        :return: Updated dataset
        """
        return dataset

    def check_prerequisites(
        self,
    ):
        """
        Check to see if any prerequisite processors are missing

        :return: None
        """

    def clean_cache(self, cache_id: int):
        """
        Function to clean the internal cache filled by base_apply

        :param cache_id: key for cache
        :return: None

        """
        del self.passed_batches[cache_id]
        del self.err_stack[cache_id]

    def base_apply(self, dataset: Dataset) -> tuple[Dataset, ErrorStack]:
        """
        Core function to act on a dataset, and return an updated dataset

        :param dataset: Input dataset
        :return: Updated dataset, and any caught errors
        """
        cache_id = threading.get_ident()

        self.passed_batches[cache_id] = {}
        self.err_stack[cache_id] = ErrorStack()

        self.latest_n_input_batches = len(dataset)
        self.latest_n_input_blocks = sum(len(x) for x in dataset)
        processor_start_time = datetime.datetime.now()
        if len(dataset) > 0:
            n_cpu = min([self.max_n_cpu, len(dataset)])

            logger.info(f"Running {self.__class__.__name__} on {n_cpu} threads")

            watchdog_queue = Queue()

            workers = []

            for _ in range(n_cpu):
                # Set up a worker thread to process database load
                worker = Thread(
                    target=self.apply_to_batch, args=(watchdog_queue, cache_id)
                )
                worker.daemon = True
                worker.start()

                workers.append(worker)

            with tqdm(total=len(dataset), position=0, leave=False) as progress:
                # Set up progress bar
                self.progress[cache_id] = progress

                # Loop over batches to add to queue
                for j, batch in enumerate(dataset):
                    watchdog_queue.put(item=(j, batch))

                # Wait for the queue to empty
                watchdog_queue.join()

                self.progress[cache_id].refresh()
                self.progress[cache_id].close()

        new_dataset = []

        for key in sorted(self.passed_batches[cache_id].keys()):
            new_dataset.append(self.passed_batches[cache_id][key])

        dataset = self.update_dataset(Dataset(new_dataset))
        err_stack = self.err_stack[cache_id]

        self.clean_cache(cache_id=cache_id)

        self.latest_n_output_batches = len(dataset)
        self.latest_n_output_blocks = sum(len(x) for x in dataset)
        self.latest_error_stack = err_stack
        processor_end_time = datetime.datetime.now()
        self.processor_run_time_minutes = (
            processor_end_time - processor_start_time
        ).total_seconds() / 60.0
        return dataset, err_stack

    def apply_to_batch(self, queue, cache_id: int):
        """
        Function to run self.apply on a batch in the queue, catch any errors, and then
        update the internal cache with the results.

        :param queue: python threading queue
        :param cache_id: key for cache
        :return: None
        """
        while True:
            j, batch = queue.get()
            try:
                batch = self.apply(batch)
                self.passed_batches[cache_id][j] = batch
            except NoncriticalProcessingError as exc:
                err = self.generate_error_report(exc, batch)
                logger.error(err.generate_log_message())
                self.err_stack[cache_id].add_report(err)
                self.passed_batches[cache_id][j] = batch
            except Exception as exc:  # pylint: disable=broad-except
                err = self.generate_error_report(exc, batch)
                logger.error(err.generate_log_message())
                self.err_stack[cache_id].add_report(err)

            self.progress[cache_id].update(1)
            self.progress[cache_id].refresh()

            queue.task_done()

    def apply(self, batch: DataBatch):
        """
        Function applying the processor to a
        :class:`~mirar.data.base_data.DataBatch`.
        Also updates the processing history.

        :param batch: input data batch
        :return: updated data batch
        """
        batch = self._apply(batch)
        batch = self._update_processing_history(batch)
        self.check_duplicates(batch)
        return batch

    def _apply(self, batch: DataBatch) -> DataBatch:
        """
        Core function to update the :class:`~mirar.data.base_data.DataBatch`

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
        :class:`~mirar.data.base_data.DataBlock` object in a
        :class:`~mirar.data.base_data.DataBatch`.

        :param batch: Input data batch
        :return: Updated data batch
        """
        for i, data_block in enumerate(batch):
            data_block[PROC_HISTORY_KEY] += self.base_key + ","
            data_block["REDUCER"] = getpass.getuser()
            data_block["REDMACH"] = socket.gethostname()
            data_block["REDTIME"] = str(datetime.datetime.now())
            data_block["REDSOFT"] = PACKAGE_NAME
            batch[i] = data_block
        return batch

    def check_duplicates(self, batch: DataBatch):
        """
        Function to check for duplicate names in a batch

        :param batch: Input data batch
        :return: None
        """
        names = [x.get_name() for x in batch]
        if len(names) != len(set(names)):
            err = f"Duplicate names in {self.__class__.__name__}"
            logger.error(err)
            raise DuplicationError(err)


class CleanupProcessor(BaseProcessor, ABC):
    """
    Processor which 'cleans up' by deleting empty batches
    """

    def update_dataset(self, dataset: Dataset) -> Dataset:
        # Remove empty dataset
        new_dataset = Dataset([x for x in dataset.get_batches() if len(x) > 0])
        return new_dataset


class ImageHandler:
    """
    Base class for handling images
    """

    @staticmethod
    def open_fits(path: str | Path) -> Image:
        """
        Opens a fits file, and returns an Image object

        :param path: Path of image
        :return: Image object
        """
        path = str(path)
        data, header = open_fits(path)
        if RAW_IMG_KEY not in header:
            header[RAW_IMG_KEY] = path
        if BASE_NAME_KEY not in header:
            header[BASE_NAME_KEY] = Path(path).name
        return Image(data=data, header=header)

    @staticmethod
    def save_fits(
        image: Image,
        path: str | Path,
        compress: bool = False,
    ):
        """
        Save an Image to path

        :param image: Image to save
        :param path: path
        :param compress: whether to compress the fits file
        :return: None
        """
        save_fits(image, path, compress=compress)

    def save_mask_image(
        self, image: Image, img_path: Path, compress: bool = False
    ) -> Path:
        """
        Saves a mask image, following the astromatic software convention of
        masked value = 0. and non-masked value = 1.

        :param image: Science image
        :param img_path: Path of parent image
        :param compress: Whether to compress the mask image
        :return: Path of mask image
        """
        mask_path = get_mask_path(img_path)
        header = image.get_header()

        mask = image.get_mask()
        if LATEST_WEIGHT_SAVE_KEY in image.header:

            path = Path(image.header[LATEST_WEIGHT_SAVE_KEY])
            if path.exists():
                weight_data = self.open_fits(
                    image.header[LATEST_WEIGHT_SAVE_KEY]
                ).get_data()
                mask = mask * weight_data
            else:
                logger.warning(
                    f"Could not find weight file {image.header[LATEST_WEIGHT_SAVE_KEY]}"
                )
        self.save_fits(Image(mask.astype(float), header), mask_path, compress=compress)

        return mask_path

    @staticmethod
    def get_hash(image_batch: ImageBatch):
        """
        Get a unique hash for an image batch

        :param image_batch: image batch
        :return: unique hash for that batch
        """
        key = "".join(
            sorted([x[BASE_NAME_KEY] + x[PROC_HISTORY_KEY] for x in image_batch])
        )
        return hashlib.sha1(key.encode()).hexdigest()


class BaseImageProcessor(BaseProcessor, ImageHandler, ABC):
    """
    Base processor handling images in/images out
    """

    def _apply(self, batch: ImageBatch) -> ImageBatch:
        return self._apply_to_images(batch)

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        raise NotImplementedError


class ProcessorWithCache(BaseImageProcessor, ABC):
    """
    Image processor with cached images associated to it, e.g a master flat
    """

    def __init__(
        self,
        try_load_cache: bool = True,
        write_to_cache: bool = True,
        overwrite: bool = True,
        cache_sub_dir: str = CAL_OUTPUT_SUB_DIR,
        cache_image_name_header_keys: str | list[str] | None = None,
    ):
        super().__init__()
        self.try_load_cache = try_load_cache
        self.write_to_cache = write_to_cache
        self.overwrite = overwrite
        self.cache_sub_dir = cache_sub_dir
        self.cache_image_name_header_keys = cache_image_name_header_keys

    def select_cache_images(self, images: ImageBatch) -> ImageBatch:
        """
        Select the appropriate cached image for the batch

        :param images: images to process
        :return: cached images to use
        """
        raise NotImplementedError

    def get_cache_path(self, images: ImageBatch) -> Path:
        """
        Gets path for saving/loading cached image

        :param images: images to process
        :return: cache path
        """

        file_name = self.get_cache_file_name(images)

        output_path = get_output_path(
            base_name=file_name, dir_root=self.cache_sub_dir, sub_dir=self.night_sub_dir
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        return output_path

    def get_cache_file_name(self, images: ImageBatch) -> str:
        """
        Get unique cache name for images

        :param images: images to process
        :return: unique hashed name
        """
        logger.debug(f"Images are {images}")
        cache_images = self.select_cache_images(images)
        cache_image_str = ""
        if (self.cache_image_name_header_keys is not None) and (len(cache_images) > 0):
            if isinstance(self.cache_image_name_header_keys, str):
                self.cache_image_name_header_keys = [self.cache_image_name_header_keys]

            cache_image_str = "_".join(
                [
                    str(cache_images[0].header[x])
                    for x in self.cache_image_name_header_keys
                ]
            )
        return f"{self.base_key}_{cache_image_str}_{self.get_hash(cache_images)}.fits"

    def get_cache_file(self, images: ImageBatch) -> Image:
        """
        Return the appropriate cached image for the batch

        :param images: images to process
        :return: cached image to use
        """

        path = self.get_cache_path(images)

        exists = path.exists()

        if np.logical_and(self.try_load_cache, exists):
            logger.debug(f"Loading cached file {path}")
            return self.open_fits(path)

        image = self.make_image(images)

        if self.write_to_cache:
            if np.sum([not exists, self.overwrite]) > 0:
                self.save_fits(image, path)

        return image

    def make_image(self, images: ImageBatch) -> Image:
        """
        Make a cached image (e.g master flat)

        :param images: images to use
        :return: cached image
        """
        raise NotImplementedError


class ProcessorPremadeCache(ProcessorWithCache, ABC):
    """
    Processor with pre-made master image
    """

    def __init__(
        self, master_image_path_generator: Callable[[ImageBatch], Path], *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.master_image_path_generator = master_image_path_generator

    def get_cache_path(self, images: ImageBatch) -> Path:
        """
        Gets path for saving/loading cached image

        :param images: Images to process
        :return: Path to cached image
        """
        return self.master_image_path_generator(images)


class BaseSourceGenerator(CleanupProcessor, ImageHandler, ABC):
    """
    Base CandidateGenerator processor (image batch in, source batch out)
    """

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def _apply(self, batch: ImageBatch) -> SourceBatch:
        source_batch = self._apply_to_images(batch)

        if len(source_batch) == 0:
            msg = "No sources found in image batch"
            logger.warning(msg)

        for batch in source_batch:
            cols = batch.get_data().columns
            for field in core_source_fields:
                if field not in cols:
                    raise MissingCoreFieldError(
                        f"Field {field} not found in source table. "
                        f"Available fields are {cols}."
                    )

        return source_batch

    def _apply_to_images(self, batch: ImageBatch) -> SourceBatch:
        raise NotImplementedError

    def get_metadata(self, image: Image) -> dict:
        """
        Get metadata from image

        :param image: Image to get metadata from
        :return: Metadata dictionary
        """
        metadata = {}

        for key in image.keys():
            if key != "COMMENT":
                metadata[key] = image[key]

        return metadata


class BaseSourceProcessor(BaseProcessor, ABC):
    """
    Base dataframe processor (Source batch in, source batch out)
    """

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses[cls.base_key] = cls

    def _apply(self, batch: SourceBatch) -> SourceBatch:
        return self._apply_to_sources(batch)

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        raise NotImplementedError

    @staticmethod
    def generate_super_dict(metadata: dict, source_row: pd.Series) -> dict:
        """
        Generate a dictionary of metadata and candidate row, with lower case keys

        :param metadata: Metadata for the source table
        :param source_row: Individual row of the source table
        :return: Combined dictionary
        """
        super_dict = {key.lower(): val for key, val in metadata.items()}
        super_dict.update(
            {key.lower(): val for key, val in source_row.to_dict().items()}
        )
        super_dict.update({key.upper(): val for key, val in super_dict.items()})
        return super_dict
