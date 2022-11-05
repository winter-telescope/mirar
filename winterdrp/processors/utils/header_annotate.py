import logging
import astropy.io.fits
import numpy as np
from winterdrp.processors.base_processor import BaseImageProcessor


logger = logging.getLogger(__name__)


class HeaderAnnotator(BaseImageProcessor):

    base_key = "header_annotator"

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

    def __str__(self) -> str:
        return f"Updates image headers by adding values for {' and '.join(self.output_key)}."

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


class HeaderEditor(BaseImageProcessor):

    base_key = "header_editor"

    def __init__(
            self,
            edit_keys: str | list[str],
            values: str | float | int | list,
    ):
        super().__init__()
        if not isinstance(edit_keys, list):
            edit_keys = [edit_keys]

        if not isinstance(values, list):
            values = [values]

        assert len(edit_keys) == len(values)
        self.edit_keys = edit_keys
        self.values = values

    def __str__(self) -> str:
        return f"Modifies image headers by updating values for {' and '.join(self.edit_keys)}."

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, header in enumerate(headers):

            for ind, key in enumerate(self.edit_keys):
                header[key] = self.values[ind]

            headers[i] = header

        return images, headers
