"""
Module to specify the input data classes for
:class:`mirar.processors.base_processor.ImageHandler`

The basic idea of the code is to pass
:class:`~mirar.data.base_data.DataBlock` objects
through a series of :class:`~wintedrp.processors.BaseProcessor` objects.
Since a given image can easily be ~10-100Mb, and there may be several hundred raw images
from a typical survey in a given night, the total data volume for these processors
could be several 10s of Gb or more. Storing these all in RAM would be very
inefficient/slow for a typical laptop or many larger processing machines.

To mitigate this, the code can be operated in **cache mode**. In that case,
after raw images are loaded, only the header data is stored in memory.
The actual image data itself is stored temporarily in as a npy file
in a dedicated cache directory, and only loaded into memory when needed.
When the data is updated, the npy file is changed.
The path of the file is a unique hash, and includes the read time of the file,
so multiple copies of an image can be read and modified independently.

In cache mode, all of the image data is temporarily stored in a cache,
and this cache can therefore reach the size of 10s of Gb.
The location of the cache is in the configurable
**output data directory**. This would increase linearly with successive code executions.
To mitigate that, and to avoid cleaning the cache by hand,
the code tries to automatically delete cache files as needed.

Python provides a default `__del__()` method for handling clean up when an object
is deleted. Images automatically delete their cache in this method. However, has a
somewhat-complicated method of 'garbage collection' (see
`the official description <https://devguide.python.org/internals/garbage-collector>`_
for more info), and it is not guaranteed that Image objects will
clean themselves.

As a fallback, when you run the code from the command line (and therefore call
__main__),  we use the standard python
`tempfile library <https://docs.python.org/3/library/tempfile.html>` to create a
temporary directory, and set this as a cache. We call the directory using `with`
context manager, ensuring that cleanup runs automatically before exiting,
even if the code crashes/raises errors. We also use `tempfile` and careful cleaning
 for the unit tests, as provided by the  base test class.
 **If you try to interact with the code in any other way, please be mindful of this
 behaviour, and ensure that you clean your cache in a responsible way!**

If you don't like this feature, you don't need to use it. Cache mode is entirely
optional, and can be disabled by setting the environment variable to false.

You can change this via an environment variable.

.. code-block:: bash

    export USE_WINTER_CACHE = false

See :doc:`usage` for more information about selecting cache mode,
and setting the output data directory.
"""
import copy
import hashlib
import logging
import threading
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io.fits import Header
from astropy.time import Time

from mirar.data.base_data import DataBatch, DataBlock
from mirar.data.cache import USE_CACHE, cache

logger = logging.getLogger(__name__)


class Image(DataBlock):
    """
    A subclass of :class:`~mirar.data.base_data.DataBlock`,
    containing an image and header.

    This class serves as input for
    :class:`~mirar.processors.base_processor.BaseImageProcessor` and
    :class:`~mirar.processors.base_processor.BaseCandidateGenerator` processors.
    """

    cache_files = []

    def __init__(self, data: np.ndarray, header: Header):
        self._data = None
        self.header = header
        super().__init__()
        if USE_CACHE:
            self.cache_path = self.get_cache_path()
            self.cache_files.append(self.cache_path)
        else:
            self.cache_path = None
        self.set_data(data=data)

    def get_cache_path(self) -> Path:
        """
        Get a unique cache path for the image (.npy file).
        This is hash, using name and time, so should be unique even
        when rerunning on the same image.

        :return: unique cache file path
        """
        base = "".join([str(Time.now()), self.get_name(), str(threading.get_ident())])
        name = f"{hashlib.sha1(base.encode()).hexdigest()}.npy"
        return cache.get_cache_dir().joinpath(name)

    def __str__(self):
        return f"<An {self.__class__.__name__} object, built from {self.get_name()}>"

    def set_data(self, data: np.ndarray):
        """
        Set the data with cache

        :param data: Updated image data
        :return: None
        """
        if USE_CACHE:
            self.set_cache_data(data)
        else:
            self.set_ram_data(data)

    def set_cache_data(self, data: np.ndarray):
        """
        Set the data with cache

        :param data: Updated image data
        :return: None
        """
        np.save(self.cache_path.as_posix(), data, allow_pickle=False)

    def set_ram_data(self, data: np.ndarray):
        """
        Set the data in RAM

        :param data: Updated image data
        :return: None
        """
        self._data = data

    def get_data(self) -> np.ndarray:
        """
        Get the image data from cache

        :return: image data (numpy array)
        """
        if USE_CACHE:
            return self.get_cache_data()

        return self.get_ram_data()

    def get_mask(self) -> np.ndarray:
        """
        Get the mask data for an image. 0 is masked, 1 is unmasked.

        :return: mask data (numpy array)
        """
        img_data = self.get_data()
        return ~np.isnan(img_data)

    def get_cache_data(self) -> np.ndarray:
        """
        Get the image data from cache

        :return: image data (numpy array)
        """
        return np.load(self.cache_path.as_posix(), allow_pickle=True)

    def get_ram_data(self) -> np.ndarray:
        """
        Get the image data from RAM

        :return: image data (numpy array)
        """
        return self._data

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

    def __del__(self):
        if self.cache_path is not None:
            self.cache_path.unlink(missing_ok=True)
            self.cache_files.remove(self.cache_path)

    def __deepcopy__(self, memo):
        new = type(self)(
            data=copy.deepcopy(self.get_data()), header=copy.deepcopy(self.get_header())
        )
        return new

    def __copy__(self):
        new = type(self)(
            data=self.get_data().__copy__(), header=self.get_header().__copy__()
        )
        return new


class ImageBatch(DataBatch):
    """
    A subclass of :class:`~mirar.data.base_data.DataBatch`,
    which contains :class:`~mirar.data.image_data.Image` objects

    To batch, de-batch, and select objects within batches, see
    :class:`~mirar.processors.utils.image_selector.ImageBatcher`,
    :class:`~mirar.processors.utils.image_selector.ImageDebatcher`, and
    :class:`~mirar.processors.utils.image_selector.ImageSelector`.
    """

    data_type = Image

    def __init__(self, batch: Optional[list[Image] | Image] = None):
        super().__init__(batch=batch)

    def append(self, item: Image):
        self._append(item)

    def get_batch(self) -> list[Image]:
        """Returns the :class:`~mirar.data.image_data.ImageBatch`
        items within the batch

        :return: list of :class:`~mirar.data.image_data.Image` objects
        """
        return self.get_data_list()
