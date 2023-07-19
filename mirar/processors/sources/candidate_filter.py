"""
Base class for source filters
"""
import logging

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
