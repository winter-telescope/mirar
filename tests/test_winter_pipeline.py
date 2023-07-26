"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.284109115600586,
    "ZP_2.0_std": 0.11124932020902634,
    "ZP_2.0_nstars": 794,
    "ZP_3.0": 23.77044677734375,
    "ZP_3.0_std": 0.09785334020853043,
    "ZP_3.0_nstars": 798,
    "ZP_4.0": 23.946128845214844,
    "ZP_4.0_std": 0.10347436368465424,
    "ZP_4.0_nstars": 800,
    "ZP_5.0": 23.99557876586914,
    "ZP_5.0_std": 0.10990961641073227,
    "ZP_5.0_nstars": 800,
    "ZP_6.0": 24.00973129272461,
    "ZP_6.0_std": 0.11433095484972,
    "ZP_6.0_nstars": 800,
    "ZP_7.0": 24.016685485839844,
    "ZP_7.0_std": 0.11775454133749008,
    "ZP_7.0_nstars": 799,
    "ZP_8.0": 24.01970672607422,
    "ZP_8.0_std": 0.12098255008459091,
    "ZP_8.0_nstars": 798,
    "ZP_AUTO": 24.017864227294922,
    "ZP_AUTO_std": 0.12125010788440704,
    "ZP_AUTO_nstars": 803,
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
