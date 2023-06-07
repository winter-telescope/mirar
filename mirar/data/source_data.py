"""
Module for SourceTable objects, and their corresponding SourceBatches
"""
from typing import Optional

import pandas as pd

from mirar.data.base_data import DataBatch, DataBlock


class SourceTable(DataBlock):
    """
    Data class for SourceTables, a type data block  based around
    sources detected in an image
    """

    def __init__(self, source_list: pd.DataFrame, metadata: dict):
        self.source_list = source_list
        self.metadata = metadata
        super().__init__()

    def get_data(self) -> pd.DataFrame:
        """
        Get the table of sources

        :return: source dataframe
        """
        return self.source_list

    def set_data(self, source_list: pd.DataFrame):
        """
        Set the table of sources

        :param source_list: new source list
        :return: None
        """
        self.source_list = source_list

    def get_metadata(self) -> dict:
        """
        Get the metadata associated with the source table

        :return: metadata
        """
        return self.metadata

    def __getitem__(self, item):
        return self.metadata.__getitem__(item)

    def __setitem__(self, key, value):
        self.metadata.__setitem__(key, value)

    def keys(self):
        """
        Return the metadata keys

        :return: keys
        """
        return self.metadata.keys()


class SourceBatch(DataBatch):
    """
    DataBatch class for holding SourceTables
    """

    data_type = SourceTable

    def __init__(self, batch: Optional[list[SourceTable] | SourceTable] = None):
        super().__init__(batch=batch)

    def append(self, item: SourceTable):
        self._append(item)

    def get_batch(self) -> list[SourceTable]:
        return self.get_data_list()
