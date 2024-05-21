"""
Module with classes to write a source table to CSV
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from mirar.data import SourceBatch
from mirar.paths import BASE_NAME_KEY, base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class CSVExporter(BaseSourceProcessor):
    """
    Class to export a source table to CSV
    """

    base_key = "CSVEXPORT"

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
        return "Processor to save sources to csv files."

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_list in batch:
            source_table = source_list.get_data()
            metadata = source_list.get_metadata()

            df_list = []
            for _, source_row in source_table.iterrows():
                super_dict = self.generate_super_dict(metadata, source_row)
                df_list.append(super_dict)

            output_dir = get_output_dir(
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )

            output_dir.mkdir(parents=True, exist_ok=True)
            csv_path = output_dir.joinpath(
                Path(metadata[BASE_NAME_KEY]).with_suffix(".csv").name
            )

            logger.debug(f"Writing source table to {csv_path}")

            df = pd.DataFrame(df_list)
            if self.export_keys is not None:
                df = df[self.export_keys]

            df.to_csv(csv_path, index=False)

        return batch
