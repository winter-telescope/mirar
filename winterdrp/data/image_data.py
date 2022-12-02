"""
Module to specify the input data classes for
:class:`winterdrp.processors.base_processor.ImageHandler`
"""
import logging

import numpy as np
from astropy.io.fits import Header

from winterdrp.data.base_data import DataBatch, DataBlock

logger = logging.getLogger(__name__)


class Image(DataBlock):
    """
    A subclass of :class:`~winterdrp.data.base_data.DataBlock`,
    containing an image and header.

    This class serves as input for
    :class:`~winterdrp.processors.base_processor.BaseImageProcessor` and
    :class:`~winterdrp.processors.base_processor.BaseCandidateGenerator` processors.
    """

    def __init__(self, data: np.ndarray, header: Header):
        self._data = None
        self.header = header
        super().__init__()
        self.set_data(data)

    def __str__(self):
        return (
            f"<An {self.__class__.__name__} object, " f"built from {self.get_name()}>"
        )

    def get_data(self) -> np.ndarray:
        """
        Get the image data

        :return: image data (numpy array)
        """
        return self._data

    def set_data(self, data: np.ndarray):
        """
        Set the data

        :param data: Updated image data
        :return: None
        """
        self._data = data

    def get_header(self) -> Header:
        """
        Get the image header

        :return: astropy Header
        """
        return self.header

    def set_header(self, header: Header):
        """
        Update the header

        :param header: updated header
        :return: None
        """
        self.header = header

    def __getitem__(self, item):
        return self.header.__getitem__(item)

    def __setitem__(self, key, value):
        self.header.__setitem__(key, value)

    def keys(self):
        """
        Get the header keys

        :return: Keys of header
        """
        return self.header.keys()


class ImageBatch(DataBatch):
    """
    A subclass of :class:`~winterdrp.data.base_data.DataBatch`,
    which contains :class:`~winterdrp.data.image_data.Image` objects
    """

    data_type = Image

    def __init__(self, batch: list[Image] | Image = None):
        super().__init__(batch=batch)

    def append(self, item: Image):
        self._append(item)

    def __str__(self):
        return (
            f"<An {self.__class__.__name__} object, "
            f"containing {[x.get_name() for x in self.get_batch()]}>"
        )

    def get_batch(self) -> list[Image]:
        """Returns the :class:`~winterdrp.data.image_data.ImageBatch`
        items within the batch

        :return: list of :class:`~winterdrp.data.image_data.Image` objects
        """
        return self.get_data_list()
