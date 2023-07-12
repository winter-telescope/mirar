"""
Module for adding metadata to Image headers
"""
import logging

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class HeaderAnnotator(BaseImageProcessor):
    """
    Processor for adding metadata to Image headers
    """

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
        return (
            f"Updates image headers by adding values "
            f"for '{self.output_key}' using '{' and '.join(self.input_keys)}'."
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for i, image in enumerate(batch):
            new_val = ""
            for key in self.input_keys:
                new_val += str(image[key])

            image[self.output_key] = new_val
            batch[i] = image

        return batch


class HeaderEditor(BaseImageProcessor):
    """
    Processor for modifying metadata in Image headers
    """

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
        return (
            f"Modifies image headers by "
            f"updating values for {' and '.join(self.edit_keys)}."
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for i, image in enumerate(batch):
            for ind, key in enumerate(self.edit_keys):
                image[key] = self.values[ind]

            batch[i] = image

        return batch
