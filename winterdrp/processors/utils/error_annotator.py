import astropy.io.fits
import numpy as np
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import raw_img_key
import logging
from winterdrp.errors import ErrorStack
from winterdrp.paths import proc_fail_key
from pathlib import Path

logger = logging.getLogger(__name__)


class ErrorStackAnnotator(BaseImageProcessor):

    base_key = "errorannotate"

    def __init__(
            self,
            errorstack: ErrorStack,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.errorstack = errorstack
        self.image_dict = self.unpack_errorstack()

    def unpack_errorstack(self) -> dict:
        image_dict = dict()

        all_reports = self.errorstack.get_all_reports()

        for error_report in all_reports:
            images = error_report.contents
            name = error_report.get_error_name()

            for image in images:
                if image not in image_dict.keys():
                    image_dict[image] = [name]
                else:
                    image_dict[image].append(name)

        return image_dict

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, header in enumerate(headers):

            base_name = str(Path(header[raw_img_key]).name)

            if base_name in self.image_dict.keys():
                header[proc_fail_key] += ",".join(self.image_dict[base_name])

            headers[i] = header

        return images, headers
