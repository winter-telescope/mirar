"""
Module to modify an image using a custom user-defined function
"""

from typing import Callable

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor


class CustomImageBatchModifier(BaseImageProcessor):
    """
    Class to modify an image using a custom user-defined function
    """

    base_key = "custom_image_modifier"

    def __init__(
        self,
        image_batch_modifier: Callable[[ImageBatch], ImageBatch],
    ):
        super().__init__()
        self.image_batch_modifier = image_batch_modifier

    def __str__(self):
        return (
            f"Processor to modify image batches using "
            f"'{self.image_batch_modifier.__name__}' function."
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        new_batch = self.image_batch_modifier(batch)
        return new_batch
