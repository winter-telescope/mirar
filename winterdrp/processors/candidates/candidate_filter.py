import logging

import pandas as pd

from winterdrp.data import SourceBatch
from winterdrp.processors.base_processor import BaseDataframeProcessor

logger = logging.getLogger(__name__)


class FilterCandidates(BaseDataframeProcessor):

    base_key = "filter"

    def __init__(self, *args, **kwargs):
        super(FilterCandidates, self).__init__(*args, **kwargs)

    def _apply_to_candidates(self, batch: SourceBatch) -> SourceBatch:
        raise NotImplementedError
