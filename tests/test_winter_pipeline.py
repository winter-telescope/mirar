"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.27437400817871,
    "ZP_2.0_std": 0.1187213808298111,
    "ZP_2.0_nstars": 757,
    "ZP_3.0": 23.757644653320312,
    "ZP_3.0_std": 0.10885193943977356,
    "ZP_3.0_nstars": 760,
    "ZP_4.0": 23.92974853515625,
    "ZP_4.0_std": 0.11107434332370758,
    "ZP_4.0_nstars": 761,
    "ZP_5.0": 23.97742462158203,
    "ZP_5.0_std": 0.11933239549398422,
    "ZP_5.0_nstars": 766,
    "ZP_6.0": 23.993120193481445,
    "ZP_6.0_std": 0.12369756400585175,
    "ZP_6.0_nstars": 766,
    "ZP_7.0": 24.002582550048828,
    "ZP_7.0_std": 0.12751182913780212,
    "ZP_7.0_nstars": 766,
    "ZP_8.0": 24.00574493408203,
    "ZP_8.0_std": 0.13095952570438385,
    "ZP_8.0_nstars": 766,
    "ZP_AUTO": 24.003620147705078,
    "ZP_AUTO_std": 0.13136684894561768,
    "ZP_AUTO_nstars": 769,
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
