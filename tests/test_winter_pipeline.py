"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.28829002380371,
    "ZP_2.0_std": 0.10393648594617844,
    # "ZP_2.0_nstars": 630,
    "ZP_3.0": 23.775022506713867,
    "ZP_3.0_std": 0.09419255703687668,
    # "ZP_3.0_nstars": 637,
    "ZP_4.0": 23.957653045654297,
    "ZP_4.0_std": 0.09448163956403732,
    # "ZP_4.0_nstars": 637,
    "ZP_5.0": 24.012550354003906,
    "ZP_5.0_std": 0.0988752618432045,
    # "ZP_5.0_nstars": 638,
    "ZP_6.0": 24.029869079589844,
    "ZP_6.0_std": 0.10236646234989166,
    # "ZP_6.0_nstars": 639,
    "ZP_7.0": 24.03977394104004,
    "ZP_7.0_std": 0.10459256917238235,
    # "ZP_7.0_nstars": 638,
    "ZP_8.0": 24.04380989074707,
    "ZP_8.0_std": 0.10918440669775009,
    # "ZP_8.0_nstars": 639,
    "ZP_AUTO": 24.041522979736328,
    "ZP_AUTO_std": 0.10777240991592407,
    # "ZP_AUTO_nstars": 640,
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
