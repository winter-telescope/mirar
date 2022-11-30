import logging
import os
from abc import ABC
from threading import Thread
from queue import Queue

from winterdrp.errors import ErrorStack, ErrorReport, NoncriticalProcessingError
from winterdrp.data.base_data import PseudoList, DataBatch, Dataset
from winterdrp.processors.base_processor import BaseDPU, BaseProcessor

logger = logging.getLogger(__name__)


class BaseChain(BaseDPU):

    def get_child_processors(self) -> list[BaseProcessor]:
        raise NotImplementedError

    def base_apply(
            self,
            dataset: Dataset
    ) -> tuple[Dataset, ErrorStack]:
        raise NotImplementedError

    def __str__(self):
        return f"<Processor Chain containing: {self.get_child_processors()}>"

    def __getitem__(self, item):
        return self.get_child_processors().__getitem__(item)

    def __setitem__(self, key, value):
        self.get_child_processors().__setitem__(key, value)

    def __add__(self, other):
        raise NotImplementedError

    def __len__(self):
        return len(self.get_child_processors())

    def __iter__(self):
        return self.get_child_processors().__iter__()


class NestedChain(BaseChain):

    def __init__(self, chains: BaseChain | list[BaseChain]):
        if not isinstance(chains, list):
            chains = [chains]
        self.chains = chains

    def get_component_chains(self):
        return self.chains

    def get_child_processors(self):
        all_processors = []
        for chain in self.get_component_chains():
            all_processors += chain.get_child_processors()
        return all_processors

    def apply(self, batch: DataBatch) -> DataBatch:
        i = 0
        for chain in self.get_component_chains():
            for processor in chain.get_child_processors():
                logger.debug(f"Applying '{processor.__class__} (Step {i + 1}/{len(self.get_child_processors())})")
                batch = processor.apply(batch)
                i += 1
        return batch

    def base_apply(
            self,
            dataset: Dataset
    ) -> tuple[Dataset, ErrorStack]:

        err_stack = ErrorStack()

        for chain in self.get_component_chains():

            dataset, new_err_stack = chain.base_apply(
                dataset
            )
            err_stack += new_err_stack

        return dataset, err_stack

    def __add__(self, other: BaseChain):
        new = self.__class__(chains=self.chains)
        if isinstance(other, self.__class__):
            new.chains += other.chains
        else:
            new.chains += [other]
        return new

    def __iadd__(self, other):
        if isinstance(other, NestedChain):
            self.chains += other.chains
        else:
            self.chains += [other]
        return self


class SingleChain(BaseChain):

    def __init__(
            self,
            processors: BaseProcessor | list[BaseProcessor],
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        if not isinstance(processors, list):
            processors = [processors]
        self.processors = processors

    def get_child_processors(self):
        return self.processors

    def apply(self, batch: DataBatch) -> DataBatch:
        for i, processor in enumerate(self.get_child_processors()):
            logger.debug(f"Applying '{processor.__class__} (Step {i + 1}/{len(self.get_child_processors())})")
            batch = processor.apply(batch)
        return batch

    def __add__(self, other: BaseChain):
        if type(other) == self.__class__:
            new = self.__class__(processors=self.processors)
            new.processors += other.get_child_processors()
        else:
            new = NestedChain([self, other])
        return new

    def base_apply(self, dataset: Dataset) -> tuple[Dataset, ErrorStack]:
        err_stack = ErrorStack()
        for processor in self.get_child_processors():
            dataset, new_err_stack = processor.base_apply(
                dataset
            )
            err_stack += new_err_stack
        return dataset, err_stack


class ParallelChain(SingleChain):

    def __init__(
            self,
            processors: BaseProcessor | list[BaseProcessor],
            n_cpu: int = max(int(os.cpu_count() / 2), 1),
            *args,
            **kwargs
    ):
        super().__init__(processors, *args, **kwargs)
        self.n_cpu = n_cpu
        self.passed_batches = self.err_stack = None

    def __str__(self):
        return f"Processor Chain to split processing onto {self.n_cpu} CPU threads, and then run: \n" \
               f" {self.get_child_processors()}"

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

        dataset = self.passed_batches
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
