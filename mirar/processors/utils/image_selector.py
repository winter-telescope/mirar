"""
Module containing processors and functions to select a subset of images from a batch
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.errors import ProcessorError
from mirar.processors.base_processor import BaseImageProcessor, CleanupProcessor

logger = logging.getLogger(__name__)


class ParsingError(KeyError, ProcessorError):
    """
    Exception arising due to errors in parsing Image headers
    """


def select_from_images(
    batch: ImageBatch,
    key: str = "target",
    target_values: str | list[str] = "science",
) -> ImageBatch:
    """
    Returns a subset of images in a batch with have values of <key> equal to
    a value in <target values>

    :param batch: image batch to sort
    :param key: header key to filter on
    :param target_values: accepted value(s) for key
    :return: image batch containing the subset of images which pass
    """

    # Enforce string in list for later matching
    if not isinstance(target_values, list):
        target_values = [str(target_values)]
    else:
        target_values = [str(x) for x in target_values]

    new_batch = ImageBatch()

    for image in batch:
        try:
            if str(image[key]) in target_values:
                new_batch.append(image)
        except KeyError as exc:
            logger.error(exc)
            raise ParsingError(exc) from exc

    return new_batch


class ImageSelector(BaseImageProcessor, CleanupProcessor):
    """
    Processor to only select a subset of images from a batch. Images can
    be selected using header keywords. For example, using:
        ImageSelector(("OBSTYPE", "SCIENCE"))
    selects Images with header["OBSTYPE"]=="SCIENCE"
    """

    base_key = "select"

    def __init__(self, *args: tuple[str, str | list[str]]):
        super().__init__()
        self.targets = args

    def __str__(self):
        reqs = []
        for target in self.targets:
            if isinstance(target[1], list):
                reqs.append(f"{target[0]} = {' or '.join(target[1])}")
            else:
                reqs.append(f"{target[0]} = {target[1]}")

        return f"Processor to select images where {', and '.join(reqs)}"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for header_key, target_values in self.targets:
            batch = select_from_images(
                batch, key=header_key, target_values=target_values
            )

        return batch


def split_images_into_batches(
    images: ImageBatch, split_key: str | list[str]
) -> Dataset:
    """
    Function to split a single :class:`~mirar.data.image_data.ImageBatch` object
    into multiple :class:`~mirar.data.base_data.DataBatch` objects.
    Each new batch will have the same value of <split_key>.
    Returns a dataset containing the new batches

    :param images: Image batch to split
    :param split_key: Key to split batch
    :return: Dataset containing new image batches
    """

    if isinstance(split_key, str):
        split_key = [split_key]

    groups = {}

    for image in images:
        uid = []

        for key in split_key:
            uid.append(str(image[key]))

        uid = "_".join(uid)

        if uid not in groups:
            groups[uid] = [image]
        else:
            groups[uid] += [image]
    logger.debug(groups)
    res = Dataset([ImageBatch(x) for x in groups.values()])

    return res


class ImageBatcher(BaseImageProcessor):
    """
    Module to split :class:`~mirar.data.image_data.ImageBatch` object
    into multiple :class:`~mirar.data.base_data.DataBatch` objects.

    Images are batched using the `split_key` argument. For example,
    you can batch by filter, like this:
        ImageBatcher(split_key="filter")
    which will return N batches for the N different filters present
    in the directory you are reducing.
    If you do not require batching at some point in your reductions,
    you can split by BASE_NAME_KEY:
        ImageBatcher(split_key=BASE_NAME_KEY)
    which returns ImageBatches of length 1, one for each file in the
    directory you're working with.
    """

    base_key = "batch"

    def __init__(self, split_key: str | list[str]):
        super().__init__()
        self.split_key = split_key

    def __str__(self) -> str:
        if isinstance(self.split_key, list):
            split = self.split_key
        else:
            split = [self.split_key]

        return (
            f"Groups images into batches, with each batch having "
            f"the same value of {' and '.join(split)}"
        )

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
    """
    Processor to group all incoming :class:`~mirar.data.image_data.ImageBatch`
    objects into a single batch.
    This is helpful if you've already batched at an earlier stage in your workflow, and
    you want to start over and batch by a different split key.
    """

    base_key = "debatch"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        return batch

    def __str__(self) -> str:
        return "Processor to combine all images into a single ImageBatch"

    def update_dataset(self, dataset: Dataset) -> Dataset:
        combo_batch = ImageBatch()

        for batch in dataset:
            combo_batch += batch

        return Dataset([combo_batch])
