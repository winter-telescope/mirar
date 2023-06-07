"""
Module with classes to write a candidate table to a pandas dataframe
"""
import logging
import os
from typing import Optional

from mirar.data import SourceBatch
from mirar.paths import base_output_dir, get_output_dir, get_output_path
from mirar.processors.base_processor import BaseDataframeProcessor

logger = logging.getLogger(__name__)


class DataframeWriter(BaseDataframeProcessor):
    """
    Class to write a candidate table to a pandas dataframe
    """

    base_key = "DFWRITE"

    def __init__(
        self, output_dir_name: Optional[str] = None, output_dir: str = base_output_dir
    ):
        super().__init__()
        self.output_dir_name = output_dir_name
        self.output_dir = output_dir
        logger.debug(f"Saving candidates to {self.output_dir_name}")

    def __str__(self) -> str:
        return "Processor to save candidates as a pickle file. "

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        try:
            os.makedirs(
                get_output_dir(
                    dir_root=self.output_dir_name,
                    sub_dir=self.night_sub_dir,
                    output_dir=self.output_dir,
                )
            )
        except OSError:
            pass

        for source_list in batch:
            candidate_table = source_list.get_data()

            df_basepath = os.path.basename(
                candidate_table.loc[0]["diffimname"]
            ).replace(".fits", ".candidates.pkl")
            df_path = get_output_path(
                df_basepath,
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )
            logger.info(f"Writing dataframe to {df_path}")
            candidate_table.to_pickle(df_path)

        return batch
