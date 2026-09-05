"""
A variant of CustomImageBatchModifier which also passes the night sub-directory to the wrapped function.
"""

import logging
from typing import Callable

from mirar.data import ImageBatch
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class NightAwareImageBatchModifier(BaseImageProcessor):
    """
    Like CustomImageBatchModifier, but the wrapped function receives the
    night sub-directory as well as the batch. Lets the function build
    and configure its own sub-processors, which need to know the night in
    order to resolve their output paths.
    """

    base_key = "night_aware_image_modifier"

    def __init__(
        self,
        image_batch_modifier: Callable[[ImageBatch, str], ImageBatch],
    ):
        super().__init__()
        self.image_batch_modifier = image_batch_modifier

    def description(self):
        return (
            f"Processor to modify image batches using "
            f"'{self.image_batch_modifier.__name__}' function."
        )

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        return self.image_batch_modifier(batch, self.night_sub_dir)
