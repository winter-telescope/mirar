import os
import pandas as pd
from winterdrp.processors.base_processor import BaseDataframeProcessor
import logging

logger = logging.getLogger(__name__)


class RegionsWriter(BaseDataframeProcessor):

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
        regions_path = os.path.join(self.output_dir_name, candidate_table.iloc[0]['diffimname'].replace('.fits', '.reg'))

        logger.info(f'Writing regions path to {regions_path}')
        with open(f"{regions_path}", 'w') as f:
            f.write('image\n')
            for ind in range(len(candidate_table)):
                row = candidate_table.iloc[ind]
                f.write(f"CIRCLE({row['X_IMAGE']+1},{row['Y_IMAGE']+1},5)\n")

        return candidate_table
