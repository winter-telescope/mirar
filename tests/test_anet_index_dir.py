"""
Module for testing the astrometry.net index directories and their coverage
"""

import logging
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

from astropy.io import fits

from mirar.processors.astrometry.anet.anet import (
    ANET_INDEX_DIR_KEY,
    RADIANS_TO_ARCMIN,
    AstrometryNetIndexError,
    check_index_coverage,
    get_field_width,
    get_index_dir,
    get_index_scales,
    get_scale_bounds,
)
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

SERIES = "5200"

ANET_LOGGER = "mirar.processors.astrometry.anet.anet"

# 1024 pixels at 0.35 arcsec/pixel is just under 6 arcmin, as for SEDMv2
IMAGE_SHAPE = (1024, 1024)
SEDMV2_SCALE_BOUNDS = [0.3, 0.35]


def write_index_file(index_dir: Path, name: str, low: float, high: float):
    """
    Write an index file holding nothing but the scale range in its header

    :param index_dir: Directory to write into
    :param name: Name of the index file
    :param low: Lowest quad scale the file holds, in arcminutes
    :param high: Highest quad scale the file holds, in arcminutes
    :return: None
    """
    header = fits.Header()
    header["SCALE_L"] = low / RADIANS_TO_ARCMIN
    header["SCALE_U"] = high / RADIANS_TO_ARCMIN
    fits.PrimaryHDU(header=header).writeto(index_dir / name)


class TestAnetScaleBounds(BaseTestCase):
    """Class for testing the scale bounds and field widths"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)

    def test_scale_bounds(self):
        """
        Test that bounds are set either side of the pixel scale

        :return: None
        """
        self.logger.info("\n\n Testing anet scale bounds \n\n")

        self.assertEqual(get_scale_bounds(0.24), [0.18, 0.30])

        for value, expected in zip(get_scale_bounds(0.2, margin=0.5), [0.1, 0.3]):
            self.assertAlmostEqual(value, expected, delta=0.0001)

    def test_field_width(self):
        """
        Test that a field width is derived from each set of scale units

        :return: None
        """
        # Arcseconds per pixel, so 1024 pixels at 0.35" is 5.97 arcmin
        self.assertAlmostEqual(
            get_field_width([0.3, 0.35], "app", IMAGE_SHAPE), 5.973, delta=0.01
        )
        self.assertEqual(get_field_width([4.0, 6.0], "amw", IMAGE_SHAPE), 6.0)
        self.assertEqual(get_field_width([0.05, 0.1], "degw", IMAGE_SHAPE), 6.0)
        self.assertIsNone(get_field_width([1.0], "parsecs", IMAGE_SHAPE))


class TestAnetIndexCoverage(BaseTestCase):
    """Class for testing ..func::check_index_coverage"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        # pylint: disable=consider-using-with
        self.base_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.base_dir.cleanup)
        self.index_dir = Path(self.base_dir.name) / SERIES
        self.index_dir.mkdir()

    def check_coverage(self, scale_bounds):
        """
        Run the coverage check against the index files written for this test

        :param scale_bounds: Scale bounds to check, in arcseconds per pixel
        :return: None
        """
        with patch.dict(os.environ, {ANET_INDEX_DIR_KEY: self.base_dir.name}):
            check_index_coverage(
                index_series=SERIES,
                scale_bounds=scale_bounds,
                scale_units="app",
                image_shape=IMAGE_SHAPE,
            )

    def test_index_scales(self):
        """
        Test that the scale range of each index file is read, in arcminutes

        :return: None
        """
        self.logger.info("\n\n Testing anet index coverage \n\n")

        write_index_file(self.index_dir, f"index-{SERIES}-43.fits", 2.0, 2.8)
        write_index_file(self.index_dir, f"index-{SERIES}1-04.fits", 2.8, 4.0)

        scales = get_index_scales(self.index_dir, SERIES)

        self.assertEqual(len(scales), 2)
        for value, expected in zip(scales, [(2.0, 2.8), (2.8, 4.0)]):
            self.assertAlmostEqual(value[0], expected[0], delta=0.001)
            self.assertAlmostEqual(value[1], expected[1], delta=0.001)

    def test_too_coarse(self):
        """
        Test that index files which cannot match the field are reported

        :return: None
        """
        write_index_file(self.index_dir, f"index-{SERIES}-43.fits", 8.0, 11.0)

        with self.assertLogs(ANET_LOGGER, level="WARNING") as logs:
            self.check_coverage(SEDMV2_SCALE_BOUNDS)

        self.assertIn("none of them can match", "".join(logs.output))

    def test_marginal(self):
        """
        Test that index files which barely fit the field are reported

        This is the real SEDMv2 case: a 6 arcmin field against a directory whose
        finest quads are 4 arcmin, which solves for one reduction but not another.

        :return: None
        """
        write_index_file(self.index_dir, f"index-{SERIES}-43.fits", 4.0, 5.6)

        with self.assertLogs(ANET_LOGGER, level="WARNING") as logs:
            self.check_coverage(SEDMV2_SCALE_BOUNDS)

        self.assertIn("marginal", "".join(logs.output))

    def test_fine_enough(self):
        """
        Test that index files which cover the field are not reported

        :return: None
        """
        write_index_file(self.index_dir, f"index-{SERIES}-43.fits", 2.0, 2.8)

        with self.assertNoLogs(ANET_LOGGER, level="WARNING"):
            self.check_coverage(SEDMV2_SCALE_BOUNDS)

    def test_no_scale_bounds(self):
        """
        Test that declaring no scale bounds at all is reported

        :return: None
        """
        write_index_file(self.index_dir, f"index-{SERIES}-43.fits", 2.0, 2.8)

        with self.assertLogs(ANET_LOGGER, level="WARNING") as logs:
            self.check_coverage(None)

        self.assertIn("No scale bounds", "".join(logs.output))


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
