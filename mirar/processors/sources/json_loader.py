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
    DATA_JSON_KEY,
    METADATA_JSON_KEY,
    SOURCE_SUFFIX,
)

logger = logging.getLogger(__name__)


def load_source_table_from_json(json_path: Path) -> SourceTable:
    """
    Function to load a source table from a json file

    :param json_path: Path to the json file
    :return: JSON file as a SourceTable
    """
    with open(json_path, "r", encoding="utf8") as json_f:
        res = json.load(json_f)
    metadata = res[METADATA_JSON_KEY]
    data = pd.read_json(res[DATA_JSON_KEY])
    return SourceTable(data, metadata)


class JsonSourceLoader(BaseSourceProcessor):
    """
    Class to write a candidate table to a pandas dataframe
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

    def __str__(self) -> str:
        return f"Processor to load sources from json files in {self.input_dir_name} . "

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
            source_table = load_source_table_from_json(source_path)
            if len(source_table) > 0:
                new_batch.append(source_table)
            else:
                logger.warning(f"Empty source table in {source_path}")

        return new_batch
