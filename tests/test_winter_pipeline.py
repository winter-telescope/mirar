"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.282760620117188,
    "ZP_2.0_std": 0.1225920245051384,
    "ZP_2.0_nstars": 721,
    "ZP_3.0": 23.771780014038086,
    "ZP_3.0_std": 0.11599168926477432,
    "ZP_3.0_nstars": 725,
    "ZP_4.0": 23.95162582397461,
    "ZP_4.0_std": 0.11737150698900223,
    "ZP_4.0_nstars": 724,
    "ZP_5.0": 24.003822326660156,
    "ZP_5.0_std": 0.12624092400074005,
    "ZP_5.0_nstars": 729,
    "ZP_6.0": 24.01872444152832,
    "ZP_6.0_std": 0.13025295734405518,
    "ZP_6.0_nstars": 729,
    "ZP_7.0": 24.025827407836914,
    "ZP_7.0_std": 0.1331106722354889,
    "ZP_7.0_nstars": 729,
    "ZP_8.0": 24.030405044555664,
    "ZP_8.0_std": 0.13536156713962555,
    "ZP_8.0_nstars": 730,
    "ZP_AUTO": 24.029008865356445,
    "ZP_AUTO_std": 0.13526619970798492,
    "ZP_AUTO_nstars": 730,
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
