"""
Tests for summer reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.375941621620427,
    "ZP_2.0_std": 0.07814918713767272,
    "ZP_3.0": 25.065196197780487,
    "ZP_3.0_std": 0.06722205837155901,
    "ZP_4.0": 25.445211601774158,
    "ZP_4.0_std": 0.06468754554779911,
    "ZP_5.0": 25.65819361030825,
    "ZP_5.0_std": 0.0665888647406662,
    "ZP_6.0": 25.77824513271701,
    "ZP_6.0_std": 0.06654911651556361,
    "ZP_7.0": 25.850486439021942,
    "ZP_7.0_std": 0.06707066267309605,
    "ZP_8.0": 25.89747143458705,
    "ZP_8.0_std": 0.06835361925053927,
    "ZP_AUTO": 25.96169715440812,
    "ZP_AUTO_std": 0.07550056418460695,
}

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
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        print("New Results SUMMER:")
        new_exp = "expected_zp = { \n"
        for header_key in header.keys():
            if header_key in expected_zp:
                new_exp += f'    "{header_key}": {header[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )
