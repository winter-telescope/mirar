"""
Base class for unit testing, with common cleanup method
"""

import json
import logging
import os
import sys
import tempfile
import unittest
from pathlib import Path

from mirar.data.cache import cache
from mirar.paths import TEMP_DIR, base_code_dir

logger = logging.getLogger(__name__)

expected_values_dir = base_code_dir / "tests" / "expected_values"

UPDATE_VALUES_KEY = "UPDATE_EXPECTED_VALUES"


def update_expected_values() -> bool:
    """
    Check whether the run should rewrite the cached values rather than assert

    :return: Whether to update
    """
    return os.getenv(UPDATE_VALUES_KEY, "false").lower() in ["true", "1", "yes"]


def get_test_platform() -> str:
    """
    Get the platform key under which expected test values are cached

    Results drift slightly between operating systems, so each one compares
    against its own set.

    :return: Platform key
    """
    if sys.platform.startswith("darwin"):
        return "macos"
    return "linux"


def get_expected_values_path(name: str) -> Path:
    """
    Get the path of the cached expected values for a test on this platform

    :param name: Name of the test, e.g. 'winter'
    :return: Path of the json file
    """
    return expected_values_dir / get_test_platform() / f"{name}.json"


def to_builtin(value):
    """
    Convert a numpy scalar to the equivalent python type, for serialising

    :param value: Value to convert
    :return: Converted value
    """
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(f"Cannot serialise a {type(value).__name__} to json")


def save_expected_values(name: str, values: dict):
    """
    Cache a set of expected values for a test on this platform

    :param name: Name of the test, e.g. 'winter'
    :param values: Values to cache
    :return: None
    """
    path = get_expected_values_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf8") as out_f:
        json.dump(values, out_f, indent=2, sort_keys=True, default=to_builtin)
        out_f.write("\n")
    logger.info(f"Wrote expected values to {path}")


def load_expected_values(name: str) -> dict:
    """
    Load the cached expected values for a test on this platform

    :param name: Name of the test, e.g. 'winter'
    :return: Cached values, or an empty dict if this platform has none
    """
    path = get_expected_values_path(name)

    if not path.exists():
        if update_expected_values():
            # Fall back to linux purely for the set of keys to collect
            path = expected_values_dir / "linux" / f"{name}.json"
            logger.warning(
                f"No values cached for {get_test_platform()} yet, writing a new set"
            )
        else:
            logger.warning(
                f"No expected values cached for '{name}' on {get_test_platform()}"
            )
            return {}

    with open(path, "r", encoding="utf8") as cache_f:
        return json.load(cache_f)


class BaseTestCase(unittest.TestCase):
    """Base TestCase object with additional cleanup"""

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        # pylint: disable=consider-using-with
        self.temp_dir = tempfile.TemporaryDirectory(dir=TEMP_DIR)
        cache.set_cache_dir(self.temp_dir.name)
        self.addCleanup(self.temp_dir.cleanup)

    def check_expected_values(self, name: str, values: dict):
        """
        Skip the test if no expected values are cached for this platform

        :param name: Name of the test, e.g. 'winter'
        :param values: Values loaded by :func:`~mirar.testing.load_expected_values`
        :return: None
        """
        if len(values) == 0:
            self.skipTest(
                f"No expected values cached for '{name}' on {get_test_platform()}. "
                f"Generate them by rerunning with {UPDATE_VALUES_KEY}=true, then "
                f"check the result before committing."
            )
