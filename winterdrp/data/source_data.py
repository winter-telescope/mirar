import pandas as pd
from winterdrp.data.base_data import Data, DataBatch


class SourceTable(Data):

    def __init__(
            self,
            source_list: pd.DataFrame,
    ):

        self._source_list = source_list

    def get_data(self) -> pd.DataFrame:
        return self._source_list

    def set_data(self, source_list: pd.DataFrame):
        self._source_list = source_list


class SourceBatch(DataBatch):

    def __add__(self, data: SourceTable):
        self.batch.append(data)

    def get_batch(self) -> list[SourceTable]:
        return self.batch
