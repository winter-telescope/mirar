"""Test suite for ..module::mirar.monitor module"""
import logging
import os
import unittest

from mirar.downloader.get_test_data import get_test_data_dir
from mirar.monitor.base_monitor import Monitor
from mirar.processors.utils.cal_hunter import CalRequirement
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

summer_cal_requirements = [
    CalRequirement(
        target_name="bias", required_field="EXPTIME", required_values=["0.0"]
    ),
    CalRequirement(
        target_name="flat", required_field="FILTERID", required_values=["r"]
    ),
]


class TestMonitor(BaseTestCase):
    """Class for testing ..module::mirar.monitor"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_tokens()

    def check_tokens(self):
        """Checks if tests should be skipped if required email doesn't exist.
        :raises unittest.SkipTest: If missing email skip.
        """
        if os.environ.get("SKIP_TEST_IF_NO_TOKEN") == "True":
            if os.environ.get("WATCHDOG_EMAIL", default="") == "":
                raise unittest.SkipTest("No Watchdog email, skipping test")

    def test_monitor(self):
        """Function to test ..class::Monitor realtime processing"""
        self.logger.info("Testing monitor")

        monitor = Monitor(
            pipeline="summer",
            cal_requirements=summer_cal_requirements,
            night="20220403",
            realtime_configurations="testlog",
            postprocess_configurations=["simrealtime"],
            final_postprocess_hours=0.01,
            midway_postprocess_hours=0.003,
            raw_dir="raw",
            base_raw_img_dir=get_test_data_dir(),
        )
        monitor.process_realtime()
        self.assertEqual(len(monitor.processed_science_images), 1)
        self.assertEqual(len(monitor.processed_cal_images), 1)
        self.assertEqual(len(monitor.failed_images), 0)
