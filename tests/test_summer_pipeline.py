"""
Tests for summer reduction
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.paths import get_output_dir
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.376594086837773,
    "ZP_2.0_std": 0.07882569714458094,
    "ZP_2.0_nstars": 30,
    "ZP_3.0": 25.077156103007,
    "ZP_3.0_std": 0.0654024719943491,
    "ZP_3.0_nstars": 30,
    "ZP_4.0": 25.449978212865194,
    "ZP_4.0_std": 0.061531141268491046,
    "ZP_4.0_nstars": 30,
    "ZP_5.0": 25.662833965174357,
    "ZP_5.0_std": 0.0636769631520071,
    "ZP_5.0_nstars": 30,
    "ZP_6.0": 25.781941561253866,
    "ZP_6.0_std": 0.06448472945729802,
    "ZP_6.0_nstars": 30,
    "ZP_7.0": 25.854275071970623,
    "ZP_7.0_std": 0.06542328480566681,
    "ZP_7.0_nstars": 30,
    "ZP_8.0": 25.90200685297648,
    "ZP_8.0_std": 0.06612588202264635,
    "ZP_8.0_nstars": 30,
    "ZP_AUTO": 25.966977791341147,
    "ZP_AUTO_std": 0.07517795091849588,
    "ZP_AUTO_nstars": 30,
}

pipeline = get_pipeline(
    instrument="summer", selected_configurations=["test"], night="20220402"
)


class TestSummerPipeline(BaseTestCase):
    """
    Module for testing summer pipeline
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
        Test summer pipeline
        Returns:

        """
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Cleanup - delete non-empty ouptut dir
        output_dir = get_output_dir(dir_root="summer/20220402")
        shutil.rmtree(output_dir)

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        print("New Results SUMMER:")
        new_exp = "expected_zp = { \n"
        for header_key in header.keys():
            if header_key in expected_zp:
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
