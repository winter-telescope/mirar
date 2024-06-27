"""
Module with classes to read a source table from a json file
"""

import json
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from mirar.data import SourceBatch, SourceTable
from mirar.paths import base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor
from mirar.processors.sources.json_exporter import (
    JSON_METADATA_KEY,
    JSON_SOURCE_KEY,
    JSON_SUFFIX,
)

logger = logging.getLogger(__name__)


def load_json_table(
    input_path: Path, metadata_key=JSON_METADATA_KEY, source_key=JSON_SOURCE_KEY
) -> SourceTable:
    """
    Function to load a source table from a json file.

    :param input_path: Path to the json file
    :param metadata_key: Metadata key to use
    :param source_key: Source key to use
    :return: SourceTable
    """
    with open(input_path, "r", encoding="utf8") as f:
        json_data = json.load(f)
        source_json = json_data[source_key]
        source_table = pd.DataFrame(source_json)
        metadata = json_data[metadata_key]

    source_table = SourceTable(source_list=source_table, metadata=metadata)

    return source_table


class JSONLoader(BaseSourceProcessor):
    """
    Class to convert a parquet file to a source table
    """

    base_key = "JSONLOAD"

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
            f"Processor to load sources from '{JSON_SUFFIX}' "
            f"json files in '{self.input_dir_name}' directory."
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

        source_tables = list(input_dir.glob(f"*{JSON_SUFFIX}"))
        for source_path in source_tables:
            source_table = load_json_table(source_path)
            if len(source_table) > 0:
                new_batch.append(source_table)
            else:
                logger.warning(f"Empty source table in {source_path}")

        return new_batch
