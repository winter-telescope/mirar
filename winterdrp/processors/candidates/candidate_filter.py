import logging

import pandas as pd

from winterdrp.processors.base_processor import BaseDataframeProcessor

logger = logging.getLogger(__name__)


class FilterCandidates(BaseDataframeProcessor):

    base_key = "filter"

    def __init__(self,
                 *args,
                 **kwargs):
        super(FilterCandidates, self).__init__(*args, **kwargs)

    def _apply_to_candidates(self, candidate_table: pd.DataFrame) -> pd.DataFrame:
        raise NotImplementedError
