"""
Base class for unit testing, with common cleanup method
"""
import tempfile
import unittest

from mirar.data.cache import cache
from mirar.paths import TEMP_DIR


class BaseTestCase(unittest.TestCase):
    """Base TestCase object with additional cleanup"""

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        # pylint: disable=consider-using-with
        self.temp_dir = tempfile.TemporaryDirectory(dir=TEMP_DIR)
        cache.set_cache_dir(self.temp_dir.name)
        self.addCleanup(self.temp_dir.cleanup)
