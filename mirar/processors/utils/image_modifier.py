"""
Module to modify an image using a custom user-defined function
"""
from typing import Callable

from mirar.data import Image, ImageBatch
from mirar.processors.base_processor import BaseImageProcessor


class CustomImageModifier(BaseImageProcessor):
    """
    Class to modify an image using a custom user-defined function
    """

    base_key = "custom_image_modifier"

    def __init__(
        self,
        image_modifier: Callable[[Image], Image],
    ):
        super().__init__()
        self.image_modifier = image_modifier

    def __str__(self):
        return (
            f"Processor to modify images using "
            f"'{self.image_modifier.__name__}' function."
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for i, image in enumerate(batch):
            batch[i] = self.image_modifier(image)
        return batch
