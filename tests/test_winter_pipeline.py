"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.29535484313965,
    "ZP_2.0_std": 0.10248779505491257,
    # "ZP_2.0_nstars": 708,
    "ZP_3.0": 23.78312873840332,
    "ZP_3.0_std": 0.09413045644760132,
    # "ZP_3.0_nstars": 714,
    "ZP_4.0": 23.969829559326172,
    "ZP_4.0_std": 0.09479842334985733,
    # "ZP_4.0_nstars": 715,
    "ZP_5.0": 24.03107261657715,
    "ZP_5.0_std": 0.09941477328538895,
    # "ZP_5.0_nstars": 717,
    "ZP_6.0": 24.05586051940918,
    "ZP_6.0_std": 0.10112165659666061,
    # "ZP_6.0_nstars": 714,
    "ZP_7.0": 24.068748474121094,
    "ZP_7.0_std": 0.1044672280550003,
    # "ZP_7.0_nstars": 713,
    "ZP_8.0": 24.074588775634766,
    "ZP_8.0_std": 0.1080445870757103,
    # "ZP_8.0_nstars": 712,
    "ZP_AUTO": 24.07155418395996,
    "ZP_AUTO_std": 0.10698756575584412,
    # "ZP_AUTO_nstars": 714,
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
