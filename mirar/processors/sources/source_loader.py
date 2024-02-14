"""
Module with classes to read a source table from a json file
"""

import logging
import pickle
from pathlib import Path
from typing import Optional

from mirar.data import SourceBatch, SourceTable
from mirar.paths import base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor
from mirar.processors.sources.source_exporter import SOURCE_SUFFIX

logger = logging.getLogger(__name__)


def load_source_table(input_path: Path) -> SourceTable:
    """
    Function to load a source table from a pickle file.

    WARNING: Usual rules with pickle apply. Only load files you trust.

    :param input_path: Path to the pickle file
    :return: SourceTable
    """
    with open(input_path, "rb") as pickle_f:
        res = pickle.load(pickle_f)
    return res


class SourceLoader(BaseSourceProcessor):
    """
    Class to write a candidate table to a pandas dataframe
    """

    base_key = "SRCLOAD"

    def __init__(
        self,
        input_dir_name: Optional[str] = None,
        input_dir: str | Path = base_output_dir,
    ):
        super().__init__()
        self.input_dir_name = input_dir_name
        self.input_dir = Path(input_dir)

    def __str__(self) -> str:
        return (
            f"Processor to load sources from '{SOURCE_SUFFIX}' "
            f"files in {self.input_dir_name} . "
        )

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        input_dir = get_output_dir(
            dir_root=self.input_dir_name,
            sub_dir=self.night_sub_dir,
            output_dir=self.input_dir,
        )
        input_dir.mkdir(parents=True, exist_ok=True)

        new_batch = SourceBatch()

        source_tables = list(input_dir.glob(f"*{SOURCE_SUFFIX}"))
        for source_path in source_tables:
            source_table = load_source_table(source_path)
            if len(source_table) > 0:
                new_batch.append(source_table)
            else:
                logger.warning(f"Empty source table in {source_path}")

        return new_batch
