"""
Central module for handling the cache, currently used only for storing image data.
"""

import logging
import os
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

USE_CACHE: bool = os.getenv("USE_WINTER_CACHE", "true") in ["true", "True", True]


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

        if not USE_CACHE:
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

    def set_cache_dir(self, cache_dir: Path | str):
        """
        Function to set the cache directory

        :param cache_dir: Cache dir to set
        :return: None
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def __str__(self):
        return f"A cache, with path {self.cache_dir}"


cache = Cache()
