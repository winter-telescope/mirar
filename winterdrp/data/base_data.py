
class Data:
    pass


class DataBatch:

    def __init__(self, batch=None):

        if batch is None:
            batch = []

        self.batch = batch

    def __add__(self, data: Data):
        raise NotImplementedError

    def get_batch(self) -> list[Data]:
        raise NotImplementedError

    def __getitem__(self, item):
        return self.batch.__getitem__(item)

    def __setitem__(self, key, value):
        return self.batch.__setitem__(key, value)

    def __len__(self):
        return self.batch.__len__()


class DataSet:

    def __init__(self, batches=None):

        if batches is None:
            batches = []

        self._batches = batches

    def __add__(self, batch: DataBatch):

        if len(self._batches) > 0:
            assert type(self._batches[0]) == type(batch)

        self._batches.__add__(batch)

    def get_batch(self) -> list[Data]:
        raise NotImplementedError

    def __getitem__(self, item):
        return self._batches.__getitem__(item)

    def __setitem__(self, key, value):
        return self._batches.__setitem__(key, value)

    def __len__(self):
        return self._batches.__len__()
