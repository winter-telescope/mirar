"""
Module for adding saved errors into image header metadata
"""
import logging

from mirar.data import ImageBatch
from mirar.errors import ErrorStack
from mirar.paths import PROC_FAIL_KEY
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class ErrorStackAnnotator(BaseImageProcessor):
    """
    Processor to update image headers with processing failurs
    """

    base_key = "errorannotate"

    def __init__(self, errorstack: ErrorStack, processed_images: list[str]):
        super().__init__()
        self.errorstack = errorstack
        self.image_dict = self.unpack_errorstack()
        self.processed_images = processed_images

    def unpack_errorstack(self) -> dict:
        """
        Convert an errorstack to an image-indexed dictionary

        :return: dictionary of errors
        """
        image_dict = {}

        all_reports = self.errorstack.get_all_reports()

        for error_report in all_reports:
            image_names = error_report.contents
            error_name = error_report.get_error_name()

            for image_name in image_names:
                if image_name not in image_dict:
                    image_dict[image_name] = [error_name]
                else:
                    image_dict[image_name].append(error_name)

        return image_dict

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image in batch:
            base_name = image.get_name()

            if base_name in self.image_dict:
                image[PROC_FAIL_KEY] += ",".join(self.image_dict[base_name])
            elif self.processed_images is not None:
                if base_name not in self.processed_images:
                    image[PROC_FAIL_KEY] += "Not Processed"

        return batch
