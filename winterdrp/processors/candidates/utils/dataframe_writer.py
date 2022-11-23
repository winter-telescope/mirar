import os
import pandas as pd
from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.paths import base_output_dir, get_output_dir, get_output_path
import logging
from winterdrp.data import SourceBatch

logger = logging.getLogger(__name__)


class DataframeWriter(BaseDataframeProcessor):

    def __init__(self,
                 output_dir_name: str = None,
                 output_dir: str = base_output_dir,
                 *args,
                 **kwargs):
        super(DataframeWriter, self).__init__(*args, **kwargs)
        self.output_dir_name = output_dir_name
        self.output_dir = output_dir
        logger.debug(f"Saving candidates to {self.output_dir_name}")

    def __str__(self) -> str:
        return f"Processor to save candidates as a pickle file. "

    def _apply_to_candidates(
            self,
            batch: SourceBatch,
    ) -> SourceBatch:

        try:
            os.makedirs(get_output_dir(
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir
            ))
        except OSError:
            pass

        for source_list in batch:

            candidate_table = source_list.get_data()

            df_basepath = os.path.basename(candidate_table.loc[0]['diffimname']).replace('.fits', '.candidates.pkl')
            df_path = get_output_path(df_basepath,
                                      dir_root=self.output_dir_name,
                                      sub_dir=self.night_sub_dir,
                                      output_dir=self.output_dir
                                      )
            logger.info(f'Writing dataframe to {df_path}')
            candidate_table.to_pickle(df_path)

        return batch
