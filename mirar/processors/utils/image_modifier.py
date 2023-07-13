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

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for i, image in enumerate(batch):
            batch[i] = self.image_modifier(image)
        return batch
