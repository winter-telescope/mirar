
class DataBlock:

    def __getitem__(self, item):
        raise NotImplementedError

    def __setitem__(self, key, value):
        raise NotImplementedError


class DataBatch:

    def __init__(self, batch: list[DataBlock] | DataBlock = None):

        if batch is None:
            batch = []
        elif isinstance(batch, DataBlock):
            batch = [batch]

        self._batch = batch

    def get_batch(self) -> list[DataBlock]:
        raise NotImplementedError

    def append(self, item):
        self._batch.append(item)

    def __getitem__(self, item):
        return self._batch.__getitem__(item)

    def __setitem__(self, key, value):
        return self._batch.__setitem__(key, value)

    def __add__(self, other):
        return DataBatch(self._batch + other.get_batch())

    def __iadd__(self, other):
        self._batch += other.get_batch()
        return self

    def __len__(self):
        return self._batch.__len__()

    def __iter__(self):
        return self._batch.__iter__()


class Dataset:

    def __init__(
            self,
            batches: list[DataBatch] | DataBatch = None
    ):

        if batches is None:
            batches = []
        elif isinstance(batches, DataBatch):
            batches = [batches]

        self._batches = batches

    def get_batches(self):
        return self._batches

    def append(self, batch: DataBatch):

        if len(self._batches) > 0:
            assert type(self._batches[0]) == type(batch)

        self._batches.append(batch)

    def __getitem__(self, item):
        return self._batches.__getitem__(item)

    def __setitem__(self, key, value):
        return self._batches.__setitem__(key, value)

    def __len__(self):
        return self._batches.__len__()

    def __add__(self, other):
        return DataBatch(self._batches + other.get_batches())

    def __iadd__(self, other):
        self._batches += other.get_batches()
        return self

    def __iter__(self):
        return self._batches.__iter__()

