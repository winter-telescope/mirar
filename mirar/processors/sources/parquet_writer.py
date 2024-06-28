"""
Module with classes to save a source table as a parquet file
"""

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from mirar.data import SourceBatch
from mirar.paths import BASE_NAME_KEY, base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


PARQUET_METADATA_KEY = "mirar_metadata"
PARQUET_SUFFIX = ".parquet"


class NpEncoder(json.JSONEncoder):
    """
    Class to encode numpy objects to json
    """

    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def export_parquet(source_table: pd.DataFrame, metadata: dict, parquet_path: Path):
    """
    Function to export a source table to parquet

    :param source_table: Table dataframe to export
    :param metadata: Metadata to add to the parquet file
    :param parquet_path: Output path
    :return: None
    """
    df = source_table

    # Parquet cannot export empty dataframes
    for col in df.columns:
        if isinstance(df[col].iloc[0], pd.DataFrame):
            if len(df[col].iloc[0]) == 0:
                df[col] = [[] for _ in df[col]]

    table = pa.Table.from_pandas(df)

    custom_metadata_bytes = pd.Series(metadata).to_json().encode("utf8")
    existing_metadata = table.schema.metadata
    merged_metadata = {
        **{PARQUET_METADATA_KEY: custom_metadata_bytes},
        **existing_metadata,
    }
    table = table.replace_schema_metadata(merged_metadata)
    pq.write_table(table, parquet_path)


class ParquetWriter(BaseSourceProcessor):
    """
    Class to save a source table as a parquet file
    """

    base_key = "PARQUETWRITE"

    def __init__(
        self,
        output_dir_name: Optional[str] = None,
        output_dir: str | Path = base_output_dir,
    ):
        super().__init__()
        self.output_dir_name = output_dir_name
        self.output_dir = Path(output_dir)

    def description(self) -> str:
        return (
            f"Processor to save sources to parquet files "
            f"with '{PARQUET_SUFFIX}' suffix."
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
            source_table = source_list.get_data().copy()
            metadata = source_list.get_metadata()
            parquet_path = output_dir.joinpath(
                Path(metadata[BASE_NAME_KEY]).with_suffix(PARQUET_SUFFIX).name
            )

            logger.debug(f"Writing source table to {parquet_path}")

            export_parquet(source_table, metadata, parquet_path)

        return batch
