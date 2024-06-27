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
            source_table = source_list.get_data()
            metadata = source_list.get_metadata()

            df_list = []
            for _, source_row in source_table.iterrows():
                super_dict = self.generate_super_dict({}, source_row)
                df_list.append(super_dict)

            df = pd.DataFrame(df_list)
            parsed = df.to_json(orient="records")
            json_data = {
                "sources": json.loads(parsed),
                "metadata": json.loads(pd.Series(metadata).to_json()),
            }

            output_dir = get_output_dir(
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )

            output_dir.mkdir(parents=True, exist_ok=True)
            json_path = output_dir.joinpath(
                Path(metadata[BASE_NAME_KEY]).with_suffix(".json").name
            )

            logger.debug(f"Writing source table to {json_path}")

            with open(json_path, "w", encoding="utf8") as f:
                json.dump(json_data, f)

        return batch
