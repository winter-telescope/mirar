"""
Base class for unit testing, with common cleanup method
"""
import unittest

from winterdrp.data.image_data import clean_cache


class BaseTestCase(unittest.TestCase):
    """Base TestCase object with additional cleanup"""

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        self.addCleanup(clean_cache)
