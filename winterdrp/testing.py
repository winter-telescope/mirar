"""
Base class for unit testing, with common cleanup method
"""
import tempfile
import unittest

from winterdrp.data.cache import cache


class BaseTestCase(unittest.TestCase):
    """Base TestCase object with additional cleanup"""

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        self.temp_dir = (
            tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        )
        cache.set_cache_dir(self.temp_dir.name)
        self.addCleanup(self.temp_dir.cleanup)
