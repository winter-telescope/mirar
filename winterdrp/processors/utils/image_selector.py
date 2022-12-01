import logging

import astropy.io.fits
import numpy as np

from winterdrp.data import DataBatch, Dataset, Image, ImageBatch
from winterdrp.errors import ProcessorError
from winterdrp.processors.base_processor import BaseImageProcessor, CleanupProcessor

logger = logging.getLogger(__name__)


class ParsingError(KeyError, ProcessorError):
    pass


def select_from_images(
    batch: ImageBatch,
    key: str = "target",
    target_values: str | list[str] = "science",
) -> ImageBatch:

    # Enforce string in list for later matching
    if not isinstance(target_values, list):
        target_values = [str(target_values)]
    else:
        target_values = [str(x) for x in target_values]

    new_batch = ImageBatch()

    for i, image in enumerate(batch):
        try:
            if str(image[key]) in target_values:
                new_batch.append(image)
        except KeyError as e:
            logger.error(e)
            raise ParsingError(e)

    return new_batch


class ImageSelector(BaseImageProcessor, CleanupProcessor):

    base_key = "select"

    def __init__(self, *args: tuple[str, str | list[str]], **kwargs):
        super().__init__(*args, **kwargs)
        self.targets = args

    def __str__(self):
        reqs = []
        for x in self.targets:
            if isinstance(x[1], list):
                reqs.append(f"{x[0]} = {' or '.join(x[1])}")
            else:
                reqs.append(f"{x[0]} = {x[1]}")

        return f"Processor to select images where {', and '.join(reqs)}"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:

        for (header_key, target_values) in self.targets:

            batch = select_from_images(
                batch, key=header_key, target_values=target_values
            )

        return batch


def split_images_into_batches(
    images: ImageBatch, split_key: str | list[str]
) -> Dataset:

    if isinstance(split_key, str):
        split_key = [split_key]

    groups = dict()

    for i, image in enumerate(images):
        uid = []

        for key in split_key:
            uid.append(str(image[key]))

        uid = "_".join(uid)

        if uid not in groups.keys():
            groups[uid] = [image]
        else:
            groups[uid] += [image]

    res = Dataset([ImageBatch(x) for x in groups.values()])

    return res


class ImageBatcher(BaseImageProcessor):

    base_key = "batch"

    def __init__(self, split_key: str | list[str], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.split_key = split_key

    def __str__(self) -> str:

        if isinstance(self.split_key, list):
            split = self.split_key
        else:
            split = [self.split_key]

        return f"Groups images into batches, with each batch having the same value of {' and '.join(split)}"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        return batch

    def update_dataset(self, dataset: Dataset) -> Dataset:

        new_dataset = Dataset()

        for batch in dataset:
            new = split_images_into_batches(batch, split_key=self.split_key)
            new_dataset += new

        return new_dataset


class ImageDebatcher(BaseImageProcessor):

    base_key = "debatch"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        return batch

    def update_dataset(self, dataset: Dataset) -> Dataset:

        combo_batch = ImageBatch()

        for batch in dataset:
            combo_batch += batch

        return Dataset([combo_batch])
