"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.25054359436035,
    "ZP_2.0_std": 0.13186413049697876,
    # "ZP_2.0_nstars": 453,
    "ZP_3.0": 23.74294090270996,
    "ZP_3.0_std": 0.12157108634710312,
    # "ZP_3.0_nstars": 450,
    "ZP_4.0": 23.923004150390625,
    "ZP_4.0_std": 0.12683318555355072,
    # "ZP_4.0_nstars": 453,
    "ZP_5.0": 23.978313446044922,
    "ZP_5.0_std": 0.13428601622581482,
    # "ZP_5.0_nstars": 458,
    "ZP_6.0": 24.00155258178711,
    "ZP_6.0_std": 0.13657590746879578,
    # "ZP_6.0_nstars": 459,
    "ZP_7.0": 24.013032913208008,
    "ZP_7.0_std": 0.13995112478733063,
    # "ZP_7.0_nstars": 461,
    "ZP_8.0": 24.016738891601562,
    "ZP_8.0_std": 0.14567917585372925,
    # "ZP_8.0_nstars": 465,
    "ZP_AUTO": 24.013835906982422,
    "ZP_AUTO_std": 0.14517906308174133,
    # "ZP_AUTO_nstars": 463,
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
