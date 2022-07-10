import os
import pandas as pd
from winterdrp.processors.base_processor import BaseDataframeProcessor
import logging

logger = logging.getLogger(__name__)


class DataframeWriter(BaseDataframeProcessor):

    def __init__(self,
                 output_dir_name: str = None,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.output_dir_name = output_dir_name

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        df_path = os.path.join(self.output_dir_name, candidate_table.iloc[0]['diffimname'].replace('.fits', '.candidates.pkl'))

        logger.info(f'Writing regions path to {df_path}')
        candidate_table.to_pickle(df_path)

        return candidate_table
