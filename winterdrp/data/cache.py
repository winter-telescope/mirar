"""
Central module for handling the cache, currently used only for storing image data.
"""

import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

logger = logging.getLogger(__name__)

USE_CACHE: bool = bool(os.getenv("USE_WINTER_CACHE", "true"))


class CacheError(Exception):
    """Error Relating to cache"""


class Cache:
    """
    A cache object for storing temporary data
    """

    cache_dir: Path | None = None

    def get_cache_dir(self) -> Path:
        """
        Returns the current cache dir

        :return: Cache dir
        """
        if np.logical_and(self.cache_dir is not None, USE_CACHE):
            return self.cache_dir

        if USE_CACHE:
            err = (
                "The code has been configured to not use a cache, "
                "but is now trying to access a cache."
            )
            logger.error(err)
            raise CacheError(err)

        err = (
            "No cache dir has been set. "
            "Please set that before trying to use the cache."
        )
        logger.error(err)
        raise CacheError(err)

    def set_cache_dir(self, cache_dir: TemporaryDirectory | Path | str):
        """
        Function to set the cache directory

        :param cache_dir: Cache dir to set
        :return: None
        """
        if isinstance(cache_dir, TemporaryDirectory):
            cache_dir = cache_dir.name
        else:
            logger.warning(
                f"Setting the cache to directory {cache_dir}, "
                f"rather than using a {TemporaryDirectory}. "
                f"Remember to clean this cache yourself after you are done!"
            )
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def __str__(self):
        return f"A cache, with path {self.cache_dir}"


cache = Cache()
