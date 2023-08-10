"""
Module for modifying a source table.
"""
import logging
from typing import Callable

import pandas as pd

from mirar.data import SourceBatch
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class CustomSourceModifier(BaseSourceProcessor):
    """
    Class to modify a source table based on a function
    """

    base_key = "custom_source_modifier"

    def __init__(self, modifier_function: Callable[[pd.DataFrame], pd.DataFrame]):
        super().__init__()
        self.modifier_function = modifier_function

    def __str__(self) -> str:
        return (
            f"Processor to modify a source dataframe using the"
            f" {self.modifier_function.__name__} function."
        )

    def _apply_to_sources(self, batch: SourceBatch) -> SourceBatch:
        modified_batch = SourceBatch()
        for source_list in batch:
            candidate_table = source_list.get_data()
            modified_table = self.modifier_function(candidate_table)
            if len(modified_table) == 0:
                logger.warning("Modified source table is empty")
            else:
                source_list.set_data(modified_table)
                modified_batch.append(source_list)
        return modified_batch
