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
        return "Processor to modify a source dataframe based on a function."

    def _apply_to_sources(self, batch: SourceBatch) -> SourceBatch:
        for source_list in batch:
            candidate_table = source_list.get_data()
            modified_table = self.modifier_function(candidate_table)
            source_list.set_data(modified_table)
        return batch
