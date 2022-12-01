import pandas as pd

from winterdrp.data.base_data import DataBatch, DataBlock


class SourceTable(DataBlock):
    def __init__(self, source_list: pd.DataFrame, metadata: dict):
        self.source_list = source_list
        self.metadata = metadata
        super().__init__()

    def get_data(self) -> pd.DataFrame:
        return self.source_list

    def set_data(self, source_list: pd.DataFrame):
        self.source_list = source_list

    def get_metadata(self) -> dict:
        return self.metadata

    def __getitem__(self, item):
        return self.metadata.__getitem__(item)

    def __setitem__(self, key, value):
        self.metadata.__setitem__(key, value)

    def keys(self):
        return self.metadata.keys()


class SourceBatch(DataBatch):

    data_type = SourceTable

    def __init__(self, batch: list[SourceTable] | SourceTable = None):
        super().__init__(batch=batch)

    def append(self, data: SourceTable):
        self._append(data)

    def get_batch(self) -> list[SourceTable]:
        return self.get_data_list()
