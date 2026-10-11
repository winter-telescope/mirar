"""
Tests for summer reduction
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.paths import get_output_dir
from mirar.pipelines import get_pipeline
from mirar.testing import (
    BaseTestCase,
    load_expected_values,
    save_expected_values,
    update_expected_values,
)

logger = logging.getLogger(__name__)

EXPECTED = load_expected_values("summer")
expected_zp = EXPECTED.get("expected_zp", {})

pipeline = get_pipeline(
    instrument="summer", selected_configurations=["test"], night="20220402"
)


class TestSummerPipeline(BaseTestCase):
    """
    Module for testing summer pipeline
    """

    def setUp(self):
        """
        Function to set up test
        Returns:

        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Test summer pipeline
        Returns:

        """
        self.check_expected_values("summer", EXPECTED)
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Cleanup - delete non-empty ouptut dir
        output_dir = get_output_dir(dir_root="summer/20220402")
        shutil.rmtree(output_dir)

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        if update_expected_values():
            save_expected_values(
                "summer", {"expected_zp": {x: header[x] for x in expected_zp}}
            )
            return

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )
