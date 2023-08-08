"""
Base class for source filters
"""
import logging
from typing import Callable

import pandas as pd

from mirar.data import SourceBatch
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class BaseSourceFilter(BaseSourceProcessor):
    """
    Base class for source filters
    """

    base_key = "filter"

    def _apply_to_sources(self, batch: SourceBatch) -> SourceBatch:
        raise NotImplementedError


class SourceFilterwithFunction(BaseSourceFilter):
    """
    Class to filter sources based on a function
    """

    base_key = "filter_using_function"

    def __init__(self, filter_function: Callable[[pd.DataFrame], pd.DataFrame]):
        super().__init__()
        self.filter_function = filter_function

    def __str__(self) -> str:
        return "Processor to filter sources based on a function."

    def _apply_to_sources(self, batch: SourceBatch) -> SourceBatch:
        for source_list in batch:
            candidate_table = source_list.get_data()
            filtered_table = self.filter_function(candidate_table)
            source_list.set_data(filtered_table)
        return batch
