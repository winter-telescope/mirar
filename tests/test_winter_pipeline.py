"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.267179489135742,
    "ZP_2.0_std": 0.10821429640054703,
    "ZP_2.0_nstars": 800,
    "ZP_3.0": 23.75172996520996,
    "ZP_3.0_std": 0.10133316367864609,
    "ZP_3.0_nstars": 803,
    "ZP_4.0": 23.929208755493164,
    "ZP_4.0_std": 0.10460660606622696,
    "ZP_4.0_nstars": 805,
    "ZP_5.0": 23.978759765625,
    "ZP_5.0_std": 0.11048795282840729,
    "ZP_5.0_nstars": 805,
    "ZP_6.0": 23.99363136291504,
    "ZP_6.0_std": 0.1148594543337822,
    "ZP_6.0_nstars": 806,
    "ZP_7.0": 24.00010108947754,
    "ZP_7.0_std": 0.11953883618116379,
    "ZP_7.0_nstars": 808,
    "ZP_8.0": 24.004154205322266,
    "ZP_8.0_std": 0.12182370573282242,
    "ZP_8.0_nstars": 805,
    "ZP_AUTO": 24.002948760986328,
    "ZP_AUTO_std": 0.12172471731901169,
    "ZP_AUTO_nstars": 808,
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
