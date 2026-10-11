"""
Module for testing the resolution of astrometry.net index directories
"""

import logging
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

from mirar.processors.astrometry.anet.anet import (
    ANET_INDEX_DIR_KEY,
    AstrometryNetIndexError,
    get_index_dir,
)
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

SERIES = "5200"


class TestAnetIndexDir(BaseTestCase):
    """Class for testing ..func::get_index_dir"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # pylint: disable=consider-using-with
        self.base_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.base_dir.cleanup)

    def test_index_dir(self):
        """
        Test that a series subdirectory holding index files is found

        :return: None
        """
        index_dir = Path(self.base_dir.name) / SERIES
        index_dir.mkdir()
        (index_dir / f"index-{SERIES}-43.fits").touch()

        with patch.dict(os.environ, {ANET_INDEX_DIR_KEY: self.base_dir.name}):
            self.assertEqual(get_index_dir(SERIES), index_dir)

    def test_no_base_dir(self):
        """
        Test that an unset base directory raises

        :return: None
        """
        with patch.dict(os.environ):
            os.environ.pop(ANET_INDEX_DIR_KEY, None)
            with self.assertRaises(AstrometryNetIndexError) as err:
                get_index_dir(SERIES)

        self.assertIn(ANET_INDEX_DIR_KEY, str(err.exception))

    def test_missing_series_dir(self):
        """
        Test that a base directory without the series subdirectory raises

        :return: None
        """
        with patch.dict(os.environ, {ANET_INDEX_DIR_KEY: self.base_dir.name}):
            with self.assertRaises(AstrometryNetIndexError) as err:
                get_index_dir(SERIES)

        message = str(err.exception)
        self.assertIn(SERIES, message)
        self.assertIn(self.base_dir.name, message)
        self.assertIn("http", message)

    def test_empty_series_dir(self):
        """
        Test that a series subdirectory with no index files raises

        :return: None
        """
        (Path(self.base_dir.name) / SERIES).mkdir()

        with patch.dict(os.environ, {ANET_INDEX_DIR_KEY: self.base_dir.name}):
            with self.assertRaises(AstrometryNetIndexError) as err:
                get_index_dir(SERIES)

        self.assertIn("index-52*.fits", str(err.exception))

    def test_wrong_series_files(self):
        """
        Test that a series subdirectory holding another series raises

        :return: None
        """
        index_dir = Path(self.base_dir.name) / SERIES
        index_dir.mkdir()
        (index_dir / "index-4202-00.fits").touch()

        with patch.dict(os.environ, {ANET_INDEX_DIR_KEY: self.base_dir.name}):
            with self.assertRaises(AstrometryNetIndexError) as err:
                get_index_dir(SERIES)

        self.assertIn("index-52*.fits", str(err.exception))
