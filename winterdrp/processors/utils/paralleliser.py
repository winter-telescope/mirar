import logging
import os
import threading
from threading import Thread
from queue import Queue
from winterdrp.paths import base_name_key
from winterdrp.processors.base_processor import BaseProcessor, CleanupProcessor
from winterdrp.data import DataBatch, Dataset
from winterdrp.errors import ErrorStack, NoncriticalProcessingError, ErrorReport

logger = logging.getLogger(__name__)


class ParallelProcessor(BaseProcessor):

    base_key = "parallel"

    def __init__(
            self,
            child_processors: list[BaseProcessor] | BaseProcessor,
            n_cpu: int = max(int(os.cpu_count() / 2), 1),
            *args,
            **kwargs
    ):
        super().__init__()
        self.n_cpu = n_cpu

        if not isinstance(child_processors, list):
            child_processors = [child_processors]

        self.child_processors = child_processors

        self.passed_batches = self.err_stack = None

    def __str__(self):
        return f"Processor to split processing onto {self.n_cpu} CPU threads, and then run: \n" \
               f" {self.child_processors}"

    def set_night(self, night_sub_dir: str | int = ""):
        super().set_night(night_sub_dir=night_sub_dir)
        for processor in self.child_processors:
            processor.set_night(night_sub_dir=night_sub_dir)

    def clean_cache(self):
        self.passed_batches = self.err_stack = None

    def base_apply(
            self,
            dataset: Dataset
    ) -> tuple[Dataset, ErrorStack]:

        self.passed_batches = Dataset()
        self.err_stack = ErrorStack()

        watchdog_queue = Queue()

        workers = []

        for i in range(self.n_cpu):
            # Set up a worker thread to process database load
            worker = Thread(target=self.apply_to_batch, args=(watchdog_queue,))
            worker.daemon = True
            worker.start()

            workers.append(worker)

        for i, batch in enumerate(dataset):
            watchdog_queue.put(item=batch)

        watchdog_queue.join()

        dataset = self.update_dataset(self.passed_batches)
        err_stack = self.err_stack

        self.clean_cache()

        return dataset, err_stack

    def apply_to_batch(
            self,
            q
    ):
        while True:
            batch = q.get()
            try:
                batch = self.apply(batch)
                self.passed_batches.append(batch)
            except NoncriticalProcessingError as e:
                err = self.generate_error_report(e, batch)
                logger.error(err.generate_log_message())
                self.err_stack.add_report(err)
                self.passed_batches.append(batch)
            except Exception as e:
                err = self.generate_error_report(e, batch)
                logger.error(err.generate_log_message())
                self.err_stack.add_report(err)
            q.task_done()

    def apply(self, batch: DataBatch) -> DataBatch:
        for processor in self.child_processors:
            batch = processor.apply(batch)
        return batch

    def generate_error_report(self, exception: Exception, batch: DataBatch):
        return ErrorReport(exception, self.__module__, [])



