"""
Module for reading/logging values from image headers
"""
import logging

from mirar.data import ImageBatch
from mirar.paths import BASE_NAME_KEY
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class HeaderReader(BaseImageProcessor):
    """
    Processor to extract data from image headers, and print it to the log
    """

    base_key = "header_reader"

    def __init__(
        self,
        *keys,
    ):
        super().__init__()
        self.keys = keys

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            msg = f"{image[BASE_NAME_KEY]} \n"
            for key in self.keys:
                msg += f"{key}: {image[key]} \n"
            logger.info(msg)

        return batch
