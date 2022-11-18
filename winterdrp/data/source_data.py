import pandas as pd
from winterdrp.data.base_data import Data, DataBatch


class SourceTable(Data):

    def __init__(
            self,
            source_list: pd.DataFrame,
    ):

        self._source_list = source_list
        self._metadata = dict()

    def get_data(self) -> pd.DataFrame:
        return self._source_list

    def set_data(self, source_list: pd.DataFrame):
        self._source_list = source_list

    def get_metadata(self) -> dict:
        return self._metadata

    def __getitem__(self, item):
        return self._metadata.__getitem__(item)

    def __setitem__(self, key, value):
        self._metadata.__setitem__(key, value)

    def keys(self):
        return self._metadata.keys()


class SourceBatch(DataBatch):

    def __add__(self, data: SourceTable):
        self._batch.append(data)

    def get_batch(self) -> list[SourceTable]:
        return self._batch
