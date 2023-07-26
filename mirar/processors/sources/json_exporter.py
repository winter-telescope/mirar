"""
Module with classes to write a candidate table to a pandas dataframe
"""
import json
import logging
from pathlib import Path
from typing import Optional

from mirar.data import SourceBatch, SourceTable
from mirar.paths import base_output_dir, get_output_dir, get_output_path
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)

SOURCE_SUFFIX = "_sources.json"
METADATA_JSON_KEY = "metadata"
DATA_JSON_KEY = "data"


def save_source_table_to_json(source_table: SourceTable, json_path: Path):
    """
    Function to save a source table to a json file

    :param source_table: SourceTable to save
    :param json_path: Path to save to
    :return: None
    """

    source_data = source_table.get_data().to_json()

    source_metadata = source_table.get_metadata()

    table_json = {
        DATA_JSON_KEY: source_data,
        METADATA_JSON_KEY: source_metadata,
    }

    with open(json_path, "w", encoding="utf8") as json_f:
        json.dump(table_json, json_f)


class JsonSourceWriter(BaseSourceProcessor):
    """
    Class to write a source table to a pair of json files
    """

    base_key = "JSONWRITE"

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
            f"Processor to save candidates to {self.output_dir_name} as a json file. "
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
            # Export sources
            source_table = source_list.get_data()
            assert len(source_table) > 0, "Candidate table is empty"

            old = Path(source_table.loc[0]["diffimname"])
            json_basepath = old.parent / f"{old.stem}{SOURCE_SUFFIX}"

            json_path = get_output_path(
                json_basepath.name,
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )

            logger.debug(f"Writing dataframe to {json_path}")

            save_source_table_to_json(source_list, json_path)

        return batch
