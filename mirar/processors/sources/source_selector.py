"""
Module containing processors and functions to select a subset of sources from a batch
"""

# pylint: disable=duplicate-code

import logging

from mirar.data import Dataset, SourceBatch
from mirar.paths import TARGET_KEY
from mirar.processors.base_processor import BaseSourceProcessor, CleanupProcessor
from mirar.processors.utils.image_selector import ParsingError

logger = logging.getLogger(__name__)


def select_from_sources(
    batch: SourceBatch,
    key: str = TARGET_KEY,
    target_values: str | list[str] = "science",
) -> SourceBatch:
    """
    Returns a subset of sources in a batch with have values of <key> equal to
    a value in <target values>

    :param batch: source batch to sort
    :param key: header key to filter on
    :param target_values: accepted value(s) for key
    :return: source batch containing the subset of sources which pass
    """

    # Enforce string in list for later matching
    if not isinstance(target_values, list):
        target_values = [str(target_values)]
    else:
        target_values = [str(x) for x in target_values]

    new_batch = SourceBatch()

    for source_table in batch:
        try:
            if str(source_table[key]) in target_values:
                new_batch.append(source_table)
        except KeyError as exc:
            logger.error(exc)
            raise ParsingError(exc) from exc

    return new_batch


class SourceSelector(BaseSourceProcessor, CleanupProcessor):
    """
    Processor to only select a subset of sources from a batch. Sources can
    be selected using header keywords. For example, using:
        SourceSelector(("OBSCLASS", "SCIENCE"))
    selects Sources with header["OBSCLASS"]=="SCIENCE"
    """

    base_key = "select"

    def __init__(self, *args: tuple[str, str | list[str]]):
        super().__init__()
        self.targets = args

    def description(self):
        reqs = []
        for target in self.targets:
            if isinstance(target[1], list):
                reqs.append(f"{target[0]} = {' or '.join(target[1])}")
            else:
                reqs.append(f"{target[0]} = {target[1]}")

        return f"Processor to select sources where {', and '.join(reqs)}"

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for header_key, target_values in self.targets:
            batch = select_from_sources(
                batch, key=header_key, target_values=target_values
            )

        return batch


def split_sources_into_batches(
    sources: SourceBatch, split_key: str | list[str]
) -> Dataset:
    """
    Function to split a single :class:`~mirar.data.source_data.SourceBatch` object
    into multiple :class:`~mirar.data.base_data.DataBatch` objects.
    Each new batch will have the same value of <split_key>.
    Returns a dataset containing the new batches

    :param sources: Source batch to split
    :param split_key: Key to split batch
    :return: Dataset containing new source batches
    """

    if isinstance(split_key, str):
        split_key = [split_key]

    groups = {}

    for source_table in sources:
        uid = []

        for key in split_key:
            uid.append(str(source_table[key]))

        uid = "_".join(uid)

        if uid not in groups:
            groups[uid] = [source_table]
        else:
            groups[uid] += [source_table]
    logger.debug(groups)
    res = Dataset([SourceBatch(x) for x in groups.values()])

    return res


class SourceBatcher(BaseSourceProcessor):
    """
    Module to split :class:`~mirar.data.source_data.SourceBatch` object
    into multiple :class:`~mirar.data.base_data.DataBatch` objects.

    Sources are batched using the `split_key` argument. For example,
    you can batch by filter, like this:
        SourceBatcher(split_key="filter")
    which will return N batches for the N different filters present
    in the directory you are reducing.
    If you do not require batching at some point in your reductions,
    you can split by BASE_NAME_KEY:
        SourceBatcher(split_key=BASE_NAME_KEY)
    which returns SourceBatches of length 1, one for each file in the
    directory you're working with.
    """

    base_key = "batch"

    def __init__(self, split_key: str | list[str]):
        super().__init__()
        self.split_key = split_key

    def description(self) -> str:
        if isinstance(self.split_key, list):
            split = self.split_key
        else:
            split = [self.split_key]

        return (
            f"Groups sources into batches, with each batch having "
            f"the same value of {' and '.join(split)}"
        )

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        return batch

    def update_dataset(self, dataset: Dataset) -> Dataset:
        new_dataset = Dataset()

        for batch in dataset:
            new = split_sources_into_batches(batch, split_key=self.split_key)
            new_dataset += new

        return new_dataset


class SourceDebatcher(BaseSourceProcessor):
    """
    Processor to group all incoming :class:`~mirar.data.source_data.SourceBatch`
    objects into a single batch.
    This is helpful if you've already batched at an earlier stage in your workflow, and
    you want to start over and batch by a different split key.
    """

    base_key = "debatch"

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        return batch

    def description(self) -> str:
        return "Processor to combine all sources into a single SourceBatch"

    def update_dataset(self, dataset: Dataset) -> Dataset:
        combo_batch = SourceBatch()

        for batch in dataset:
            combo_batch += batch

        return Dataset([combo_batch])


class SourceRebatcher(SourceBatcher):
    """
    Processor to regroup all incoming :class:`~mirar.data.source_data.SourceBatch`
    objects into a single batch, and then split by new keys.
    This is helpful if you've already batched at an earlier stage in your workflow, and
    you want to start over and batch by a different split key.
    """

    base_key = "rebatch"

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        return batch

    def description(self) -> str:
        if isinstance(self.split_key, list):
            split = self.split_key
        else:
            split = [self.split_key]

        return f"Processor to regroup sources into batches by {' and '.join(split)}"

    def update_dataset(self, dataset: Dataset) -> Dataset:
        combo_batch = SourceBatch()

        for batch in dataset:
            combo_batch += batch

        dataset = split_sources_into_batches(combo_batch, split_key=self.split_key)

        return dataset
