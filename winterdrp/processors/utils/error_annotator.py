import logging
from pathlib import Path

import astropy.io.fits
import numpy as np

from winterdrp.data import ImageBatch
from winterdrp.errors import ErrorStack
from winterdrp.paths import proc_fail_key, raw_img_key
from winterdrp.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class ErrorStackAnnotator(BaseImageProcessor):

    base_key = "errorannotate"

    def __init__(
        self, errorstack: ErrorStack, processed_images: list[str], *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.errorstack = errorstack
        self.image_dict = self.unpack_errorstack()
        self.processed_images = processed_images

    def unpack_errorstack(self) -> dict:
        image_dict = dict()

        all_reports = self.errorstack.get_all_reports()

        for error_report in all_reports:
            image_names = error_report.contents
            error_name = error_report.get_error_name()

            for image_name in image_names:
                if image_name not in image_dict.keys():
                    image_dict[image_name] = [error_name]
                else:
                    image_dict[image_name].append(error_name)

        return image_dict

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:

        for i, image in enumerate(batch):

            base_name = str(Path(image[raw_img_key]).name)

            if base_name in self.image_dict.keys():
                image[proc_fail_key] += ",".join(self.image_dict[base_name])
            elif self.processed_images is not None:
                if base_name not in self.processed_images:
                    image[proc_fail_key] += "Not Processed"

        return batch
