import numpy as np
from astropy.io.fits import Header
from winterdrp.data.base_data import Data, DataBatch
from winterdrp.paths import raw_img_key, base_name_key


class Image(Data):

    def __init__(
            self,
            image: np.ndarray,
            header: Header
    ):

        self._data = image
        self._header = header
        self.raw_img_path = self[raw_img_key]
        self.base_name = self[base_name_key]

    def get_data(self) -> np.ndarray:
        return self._data

    def set_data(self, data: np.ndarray):
        self._data = data

    def get_header(self) -> Header:
        return self._header

    def __getitem__(self, item):
        return self._header.__getitem__(item)

    def __setitem__(self, key, value):
        self._header.__setitem__(key, value)


class ImageBatch(DataBatch):

    def append(self, data: Image):
        self._batch.append(data)

    def get_batch(self) -> list[Image]:
        return self._batch
