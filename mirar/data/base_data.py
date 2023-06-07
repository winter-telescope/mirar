"""
This contains the base data classes for the :module:`wintedrp.processors`.

The smallest unit is a :class:`~mirar.data.base_data.DataBlock` object,
corresponding to a single image.
These :class:`~mirar.data.base_data.DataBlock` objects are grouped into
:class:`~mirar.data.base_data.DataBatch` objects.
Each :class:`~wintedrp.processors.BaseProcessor` will operate on a individual
:class:`~mirar.data.base_data.DataBatch` object.

The :class:`~mirar.data.base_data.DataBatch` objects are stored within a larger
:class:`~mirar.data.base_data.DataSet` object.
A :class:`~wintedrp.processors.BaseProcessor` will iterate over each
:class:`~mirar.data.base_data.DataBatch` in a
:class:`~mirar.data.base_data.Dataset`.
"""
import logging
from pathlib import Path
from typing import Optional, Type

from mirar.paths import BASE_NAME_KEY, RAW_IMG_KEY

logger = logging.getLogger(__name__)


class DataBlock:
    """Base unit for processing, corresponding to a single image."""

    def __init__(self):
        self.raw_img_list = [Path(x) for x in self[RAW_IMG_KEY].split(",")]
        self.base_name = self[BASE_NAME_KEY]

    def __getitem__(self, item):
        raise NotImplementedError

    def __setitem__(self, key, value):
        raise NotImplementedError

    def get_name(self) -> str:
        """Function to retrieve the :variable:`mirar.paths.BASE_NAME_KEY`
        of the parent image

        :return: Base name of parent image
        """
        return self.base_name

    def get_raw_img_list(self) -> list[Path]:
        """Function to retrieve the paths of all raw images from
        which this object is derived.
        Because of stacking, this list may include multiple entries.

        :return: List of path strings
        """
        return self.raw_img_list


class PseudoList:
    """
    Base Class for a list-like object which contains a list of data.
    Other classes inherit from this object.

    The basic idea is that this class holds all the functions
    for safely creating an object with a specified data type.

    This class also contains the relevant magic functions so that `len(x)`, `x[i] = N`,
    and `for y in x` work as intended.
    """

    @property
    def data_type(self):
        """
        Each list should take one specific data type.
        This is where that type is defined.
        """
        raise NotImplementedError()

    def __init__(self, data_list=None):
        self._datalist = []

        if data_list is None:
            data_list = []
        elif isinstance(data_list, self.data_type):
            data_list = [data_list]

        if not isinstance(data_list, list):
            err = f"Found {data_list} of type {type(data_list)}. Expected a list."
            logger.error(err)
            raise ValueError(err)

        for item in data_list:
            self.append(item)

    def get_data_list(self):
        """
        Retrieve the data list

        :return: The saved list of objects
        """
        return self._datalist

    def append(self, item):
        """
        Function to append, list-style, new objects.

        :param item: Object to be added
        :return: None
        """
        self._append(item)

    def _append(self, item):
        """
        Protected method to append, list-style, new objects.
        This function also checks the data type to ensure they are correct.

        :param item: Object to be added
        :return: None
        """

        if not isinstance(item, self.data_type):
            err = (
                f"Error appending item {item} of type {type(item)}. "
                f"Expected a {self.data_type} item"
            )
            logger.error(err)
            raise ValueError(err)

        if len(self._datalist) > 0:
            if not isinstance(item, type(self._datalist[0])):
                err = (
                    f"Error appending item {item} of type {type(item)}. "
                    f"This {self.__class__.__name__} object already contains "
                    f"data of type {type(self._datalist[0])}. "
                    f"Please ensure all data is of the same type."
                )
                logger.error(err)
                raise ValueError(err)
        self._datalist.append(item)

    def __getitem__(self, item):
        return self._datalist.__getitem__(item)

    def __setitem__(self, key, value):
        self._datalist.__setitem__(key, value)

    def __add__(self, other):
        new = self.__class__()
        for item in self.get_data_list():
            new.append(item)
        for item in other.get_data_list():
            new.append(item)
        return new

    def __iadd__(self, other):
        for item in other.get_data_list():
            self._datalist.append(item)
        return self

    def __len__(self):
        return self._datalist.__len__()

    def __iter__(self):
        return self._datalist.__iter__()


class DataBatch(PseudoList):
    """
    Base class for a collection of individual
    :class:`~mirar.data.base_data.DataBlock` objects.
    Each :class:`~mirar.data.base_data.DataBatch` will be operated on
    by a :class:`~wintedrp.processors.BaseProcessor`
    """

    @property
    def data_type(self) -> Type[DataBlock]:
        raise NotImplementedError()

    def __init__(self, batch: Optional[list[DataBlock] | DataBlock] = None):
        super().__init__(data_list=batch)

    def get_batch(self) -> list[DataBlock]:
        """Returns the :class:`~mirar.data.base_data.DataBlock`
        items within the batch

        :return: list of :class:`~mirar.data.base_data.DataBlock` objects
        """
        return self.get_data_list()

    def get_raw_image_names(self) -> list[Path]:
        """Returns the name of each parent raw image

        :return: list of raw image names
        """
        img_list = []
        for data_block in self.get_batch():
            img_list += [Path(x).name for x in data_block.get_raw_img_list()]
        return img_list

    def __str__(self):
        return (
            f"<An {self.__class__.__name__} object, "
            f"containing {[x.get_name() for x in self.get_batch()]}>"
        )


class Dataset(PseudoList):
    """
    Base class for a collection of individual
    :class:`~mirar.data.base_data.DataBatch` objects.
    A :class:`~wintedrp.processors.BaseProcessor` will iterate over these.
    """

    data_type = DataBatch

    def get_batches(self):
        """Returns the :class:`~mirar.data.base_data.DataBatch`
        items within the batch

        :return: list of :class:`~mirar.data.base_data.DataBatch` objects
        """
        return self.get_data_list()

    def __init__(self, batches: Optional[list[DataBatch] | DataBatch] = None):
        super().__init__(data_list=batches)

    def append(self, item: DataBatch):
        """
        Function to append, list-style, new objects.

        :param item: Object to be added
        :return: None
        """
        super()._append(item)
