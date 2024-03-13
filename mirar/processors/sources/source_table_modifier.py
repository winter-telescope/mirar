"""
Module for modifying a source table.
"""

import logging
from typing import Callable

from mirar.data import SourceBatch
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class CustomSourceTableModifier(BaseSourceProcessor):
    """
    Class to modify a source table based on a function
    """

    base_key = "custom_source_modifier"

    def __init__(self, modifier_function: Callable[[SourceBatch], SourceBatch]):
        super().__init__()
        self.modifier_function = modifier_function

    def __str__(self) -> str:
        return (
            f"Processor to modify a source dataframe using the"
            f" {self.modifier_function.__name__} function."
        )

    def _apply_to_sources(self, batch: SourceBatch) -> SourceBatch:
        modified_batch = self.modifier_function(batch)
        return modified_batch
