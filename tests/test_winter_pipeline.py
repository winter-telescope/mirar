"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.24860191345215,
    "ZP_2.0_std": 0.13277117908000946,
    "ZP_2.0_nstars": 470,
    "ZP_3.0": 23.736268997192383,
    "ZP_3.0_std": 0.12201324105262756,
    "ZP_3.0_nstars": 468,
    "ZP_4.0": 23.91571807861328,
    "ZP_4.0_std": 0.136481374502182,
    "ZP_4.0_nstars": 477,
    "ZP_5.0": 23.973907470703125,
    "ZP_5.0_std": 0.14150625467300415,
    "ZP_5.0_nstars": 480,
    "ZP_6.0": 23.996763229370117,
    "ZP_6.0_std": 0.14331680536270142,
    "ZP_6.0_nstars": 481,
    "ZP_7.0": 24.011491775512695,
    "ZP_7.0_std": 0.14509068429470062,
    "ZP_7.0_nstars": 481,
    "ZP_8.0": 24.017929077148438,
    "ZP_8.0_std": 0.14915655553340912,
    "ZP_8.0_nstars": 482,
    "ZP_AUTO": 24.015439987182617,
    "ZP_AUTO_std": 0.1488943099975586,
    "ZP_AUTO_nstars": 481,
}


pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230710"
)

logging.basicConfig(level=logging.DEBUG)


class TestWinterPipeline(BaseTestCase):
    """
    Module for testing winter pipeline
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
        Test winter pipeline
        Returns:

        """
        self.logger.info("\n\n Testing winter pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        # # Uncomment to print new expected ZP dict
        print("New Results:")
        new_exp = "expected_zp = { \n"
        for header_key in header.keys():
            if "ZP_" in header_key:
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
