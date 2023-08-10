"""
Module with classes to write a source table
"""
import logging
import pickle
from pathlib import Path
from typing import Optional

from mirar.data import SourceBatch, SourceTable
from mirar.paths import DIFF_IMG_KEY, base_output_dir, get_output_dir, get_output_path
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)

SOURCE_SUFFIX = "_sources.pkl"


def save_source_table(source_table: SourceTable, out_path: Path):
    """
    Function to save a source table to a json file.

    :param source_table: SourceTable to save
    :param out_path: Path to save to
    :return: None
    """
    with open(out_path, "wb") as pickle_f:
        pickle.dump(source_table, pickle_f, protocol=pickle.HIGHEST_PROTOCOL)


class SourceWriter(BaseSourceProcessor):
    """
    Class to write a source table to a pair of json files
    """

    base_key = "SRCWRITE"

    def __init__(
        self,
        output_dir_name: Optional[str] = None,
        output_dir: str | Path = base_output_dir,
    ):
        super().__init__()
        self.output_dir_name = output_dir_name
        self.output_dir = Path(output_dir)

    def __str__(self) -> str:
        return (
            f"Processor to save candidates to {self.output_dir_name} "
            f"as '{SOURCE_SUFFIX}' files."
        )

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        output_dir = get_output_dir(
            dir_root=self.output_dir_name,
            sub_dir=self.night_sub_dir,
            output_dir=self.output_dir,
        )
        output_dir.mkdir(parents=True, exist_ok=True)

        for source_list in batch:
            source_table = source_list.get_data()
            assert len(source_table) > 0, "Source table is empty"

            old = Path(source_table.loc[0][DIFF_IMG_KEY])
            base_path = old.parent / f"{old.stem}{SOURCE_SUFFIX}"

            pkl_path = get_output_path(
                base_path.name,
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )

            logger.debug(f"Writing SourceTable to {pkl_path}")

            save_source_table(source_list, pkl_path)

        return batch
