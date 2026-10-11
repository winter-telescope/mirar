"""
Module for testing the CalHunter archival calibration search
"""

import logging
from pathlib import Path

from mirar.data import ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.pipelines.summer.load_summer_image import load_raw_summer_image
from mirar.processors.utils.cal_hunter import (
    CalHunter,
    CalRequirement,
    find_required_cals,
)
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

NIGHT = "20220403"

EXPECTED_N_CALS = 10


def get_requirements() -> list[CalRequirement]:
    """
    Get the calibration requirements to search for

    :return: List of calibration requirements
    """
    return [
        CalRequirement(
            target_name="bias", required_field="EXPTIME", required_values=["0.0"]
        ),
        CalRequirement(
            target_name="flat", required_field="FILTERID", required_values=["r"]
        ),
    ]


class TestCalHunter(BaseTestCase):
    """Class for testing ..module::mirar.processors.utils.cal_hunter"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_find_required_cals(self):
        """
        Test that archival calibrations are found in a preceding night

        :return: None
        """
        self.logger.info("\n\n Testing find_required_cals \n\n")

        latest_dir = Path(get_test_data_dir()) / "summer" / NIGHT / "raw"

        new_batch = find_required_cals(
            latest_dir=latest_dir,
            night=NIGHT,
            requirements=get_requirements(),
            open_f=load_raw_summer_image,
            skip_latest_night=True,
        )

        self.assertEqual(len(new_batch), EXPECTED_N_CALS)

    def test_cal_hunter_processor(self):
        """
        Test that the CalHunter processor finds archival calibrations

        :return: None
        """
        self.logger.info("\n\n Testing CalHunter \n\n")

        processor = CalHunter(
            load_image=load_raw_summer_image,
            requirements=get_requirements(),
            input_img_dir=Path(get_test_data_dir()),
        )
        processor.set_night(night_sub_dir=f"summer/{NIGHT}")

        new_batch = processor.apply(ImageBatch())

        self.assertEqual(len(new_batch), EXPECTED_N_CALS)
