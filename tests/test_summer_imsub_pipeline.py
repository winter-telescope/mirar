"""
Module to test summer image subtraction pipeline
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.paths import get_output_dir
from mirar.pipelines.summer.summer_pipeline import SummerPipeline
from mirar.testing import (
    BaseTestCase,
    load_expected_values,
    save_expected_values,
    update_expected_values,
)

logger = logging.getLogger(__name__)

EXPECTED = load_expected_values("summer_imsub")
expected_values = EXPECTED.get("expected_values", {})


pipeline = SummerPipeline(night="20220815", selected_configurations=["test_imsub"])


class TestSummerImsubPipeline(BaseTestCase):
    """
    Class to test summer imsub pipeline
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Function to test summer image sub pipeline
        Returns:

        """
        self.check_expected_values("summer_imsub", EXPECTED)
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res), 1)

        # Cleanup
        output_dir = get_output_dir("summer/20220815")
        shutil.rmtree(output_dir)

        header = res[0][0].get_header()

        if update_expected_values():
            save_expected_values(
                "summer_imsub",
                {"expected_values": {x: header[x] for x in expected_values}},
            )
            return

        for key, value in expected_values.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )
