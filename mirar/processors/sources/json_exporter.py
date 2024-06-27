"""
Module with classes to write a source table to Json
"""

import json
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from mirar.data import SourceBatch
from mirar.paths import BASE_NAME_KEY, base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)

JSON_METADATA_KEY = "metadata"
JSON_SOURCE_KEY = "sources"
JSON_SUFFIX = ".json"


class JSONExporter(BaseSourceProcessor):
    """
    Class to export a source table to JSON
    """

    base_key = "JSONEXPORT"

    def __init__(
        self,
        output_dir_name: Optional[str] = None,
        output_dir: str | Path = base_output_dir,
        export_keys: Optional[list[str]] = None,
    ):
        super().__init__()
        self.output_dir_name = output_dir_name
        self.output_dir = Path(output_dir)
        self.export_keys = export_keys

    def description(self) -> str:
        return "Processor to save sources to json files."

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_list in batch:
            source_table = source_list.get_data().copy()
            metadata = source_list.get_metadata()

            binary_cols = [
                x
                for x in source_table.columns
                if pd.api.types.infer_dtype(source_table[x]) == "bytes"
            ]

            if len(binary_cols) > 0:
                logger.debug(
                    f"The following columns contain binary data: {binary_cols}. "
                    f"These will not be exported to json."
                )
                source_table = source_table.drop(columns=binary_cols)

            # Use pandas to convert the source table to a json object and back
            # This ensures json-able data types are used
            json_data = {
                JSON_SOURCE_KEY: json.loads(source_table.to_json(orient="records")),
                JSON_METADATA_KEY: json.loads(pd.Series(metadata).to_json()),
            }

            output_dir = get_output_dir(
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )

            output_dir.mkdir(parents=True, exist_ok=True)
            json_path = output_dir.joinpath(
                Path(metadata[BASE_NAME_KEY]).with_suffix(JSON_SUFFIX).name
            )

            logger.debug(f"Writing source table to {json_path}")

            with open(json_path, "w", encoding="utf8") as f:
                json.dump(json_data, f)

        return batch
