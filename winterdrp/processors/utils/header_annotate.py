import logging
import astropy.io.fits
import numpy as np
from winterdrp.processors.base_processor import BaseProcessor


logger = logging.getLogger(__name__)


class HeaderAnnotator(BaseProcessor):

    base_key = "header_reader"

    def __init__(
            self,
            input_keys: str | list[str],
            output_key: str,
    ):
        super().__init__()
        if not isinstance(input_keys, list):
            input_keys = [input_keys]

        self.input_keys = input_keys
        self.output_key = output_key

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, header in enumerate(headers):

            new_val = ""
            for key in self.input_keys:
                new_val += str(header[key])

            header[self.output_key] = new_val
            headers[i] = header

        return images, headers
