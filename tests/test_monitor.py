#!/usr/bin/env python
import unittest
import logging
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.monitor.base_monitor import Monitor
from winterdrp.processors.utils.cal_hunter import CalRequirement

logger = logging.getLogger(__name__)

logging.getLogger("winterdrp").setLevel("DEBUG")

summer_cal_requirements = [
    CalRequirement(target_name="bias", required_field="EXPTIME", required_values=["0.0"]),
    CalRequirement(target_name="flat", required_field="FILTERID", required_values=["r"]),
]


class TestErrors(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
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
            base_raw_img_dir=get_test_data_dir()
        )
        monitor.process_realtime()
        self.assertEqual(len(monitor.processed_science_images), 1)


