"""Test suite for ..module::mirar.monitor module"""

import logging
import os

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
        """If required tokens do not exist raise an error.
        :raises RuntimeError: If missing token.
        """
        if os.getenv("WATCHDOG_EMAIL", default="") == "":
            raise RuntimeError(
                "No email sender. Set environment variable WATCHDOG_EMAIL to test."
            )
        if os.getenv("WATCHDOG_EMAIL_PASSWORD", default="") == "":
            raise RuntimeError(
                "No email password. Set environment variable WATCHDOG_EMAIL_PASSWORD to test."
            )
        if os.getenv("WATCHDOG_EMAIL_RECIPIENTS", default="") == "":
            raise RuntimeError(
                "No email recipients. Set environment variable WATCHDOG_EMAIL_RECIPIENTS to test."
            )

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
