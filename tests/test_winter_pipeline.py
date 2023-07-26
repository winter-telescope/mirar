"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.251432418823242,
    "ZP_2.0_std": 0.13307902216911316,
    "ZP_2.0_nstars": 449,
    "ZP_3.0": 23.74408721923828,
    "ZP_3.0_std": 0.11848931759595871,
    "ZP_3.0_nstars": 443,
    "ZP_4.0": 23.923213958740234,
    "ZP_4.0_std": 0.12503516674041748,
    "ZP_4.0_nstars": 447,
    "ZP_5.0": 23.97820281982422,
    "ZP_5.0_std": 0.13318082690238953,
    "ZP_5.0_nstars": 452,
    "ZP_6.0": 24.001344680786133,
    "ZP_6.0_std": 0.1352919191122055,
    "ZP_6.0_nstars": 453,
    "ZP_7.0": 24.010135650634766,
    "ZP_7.0_std": 0.1440536230802536,
    "ZP_7.0_nstars": 460,
    "ZP_8.0": 24.016393661499023,
    "ZP_8.0_std": 0.14371676743030548,
    "ZP_8.0_nstars": 459,
    "ZP_AUTO": 24.013452529907227,
    "ZP_AUTO_std": 0.1454823613166809,
    "ZP_AUTO_nstars": 459,
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

        # Expect two datasets, for two different sub-boards
        self.assertEqual(len(res[0]), 1)
        self.assertEqual(len(res[1]), 1)

        headers = [x[0].get_header() for x in res]
        header = [x for x in headers if x["SUBCOORD"] == "0_1"][0]

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
