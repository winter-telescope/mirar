"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 22.904813766479492,
    "ZP_2.0_std": 0.12673653662204742,
    "ZP_2.0_nstars": 103,
    "ZP_3.0": 23.39082908630371,
    "ZP_3.0_std": 0.11679290980100632,
    "ZP_3.0_nstars": 103,
    "ZP_4.0": 23.578359603881836,
    "ZP_4.0_std": 0.1077757328748703,
    "ZP_4.0_nstars": 101,
    "ZP_5.0": 23.638534545898438,
    "ZP_5.0_std": 0.11371082067489624,
    "ZP_5.0_nstars": 102,
    "ZP_6.0": 23.664098739624023,
    "ZP_6.0_std": 0.11352476477622986,
    "ZP_6.0_nstars": 102,
    "ZP_7.0": 23.679277420043945,
    "ZP_7.0_std": 0.11231032013893127,
    "ZP_7.0_nstars": 102,
    "ZP_8.0": 23.688657760620117,
    "ZP_8.0_std": 0.11235523223876953,
    "ZP_8.0_nstars": 102,
    "ZP_AUTO": 23.683950424194336,
    "ZP_AUTO_std": 0.1190515086054802,
    "ZP_AUTO_nstars": 102,
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
