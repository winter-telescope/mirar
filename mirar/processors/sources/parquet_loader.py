"""
Module with classes to read a source table from a parquet file
"""

import json
import logging
from pathlib import Path
from typing import Optional

import pyarrow.parquet as pq

from mirar.data import SourceBatch, SourceTable
from mirar.paths import base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor
from mirar.processors.sources.parquet_writer import PARQUET_METADATA_KEY, PARQUET_SUFFIX

logger = logging.getLogger(__name__)


def load_parquet_table(
    input_path: Path, metadata_key=PARQUET_METADATA_KEY
) -> SourceTable:
    """
    Function to load a source table from a parquet file.

    :param input_path: Path to the parquet file
    :param metadata_key: Metadata key to use
    :return: SourceTable
    """

    tab = pq.read_table(input_path)
    metadata = json.loads(tab.schema.metadata[metadata_key.encode("utf8")])
    source_table = SourceTable(source_list=tab.to_pandas(), metadata=metadata)

    return source_table


class ParquetLoader(BaseSourceProcessor):
    """
    Class to convert a parquet file to a source table
    """

    base_key = "PRQULOAD"

    def __init__(
        self,
        input_dir_name: Optional[str] = None,
        input_dir: str | Path = base_output_dir,
    ):
        super().__init__()
        self.input_dir_name = input_dir_name
        self.input_dir = Path(input_dir)

    def description(self) -> str:
        return (
            f"Processor to load sources from '{PARQUET_SUFFIX}' "
            f"parquet files in '{self.input_dir_name}' directory."
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

        source_tables = list(input_dir.glob(f"*{PARQUET_SUFFIX}"))
        for source_path in source_tables:
            source_table = load_parquet_table(source_path)
            if len(source_table) > 0:
                new_batch.append(source_table)
            else:
                logger.warning(f"Empty source table in {source_path}")

        return new_batch
