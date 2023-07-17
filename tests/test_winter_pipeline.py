"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.12510108947754,
    "ZP_2.0_std": 0.14948974549770355,
    "ZP_2.0_nstars": 142,
    "ZP_3.0": 23.61602783203125,
    "ZP_3.0_std": 0.1341128647327423,
    "ZP_3.0_nstars": 142,
    "ZP_4.0": 23.80809211730957,
    "ZP_4.0_std": 0.14140582084655762,
    "ZP_4.0_nstars": 142,
    "ZP_5.0": 23.880441665649414,
    "ZP_5.0_std": 0.14471492171287537,
    "ZP_5.0_nstars": 142,
    "ZP_6.0": 23.910764694213867,
    "ZP_6.0_std": 0.14723333716392517,
    "ZP_6.0_nstars": 142,
    "ZP_7.0": 23.92569351196289,
    "ZP_7.0_std": 0.15130828320980072,
    "ZP_7.0_nstars": 142,
    "ZP_8.0": 23.933137893676758,
    "ZP_8.0_std": 0.1561094969511032,
    "ZP_8.0_nstars": 142,
    "ZP_AUTO": 23.93084716796875,
    "ZP_AUTO_std": 0.15383706986904144,
    "ZP_AUTO_nstars": 142,
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
