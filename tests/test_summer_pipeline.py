"""
Tests for summer reduction
"""
import logging

from winterdrp.data import Dataset, ImageBatch
from winterdrp.pipelines import get_pipeline
from winterdrp.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.379148900858567,
    "ZP_2.0_std": 0.07320184249312667,
    "ZP_2.0_nstars": 30,
    "ZP_3.0": 25.072568754704793,
    "ZP_3.0_std": 0.06529106709538676,
    "ZP_3.0_nstars": 30,
    "ZP_4.0": 25.449809682718914,
    "ZP_4.0_std": 0.06168474184669056,
    "ZP_4.0_nstars": 30,
    "ZP_5.0": 25.66054034665426,
    "ZP_5.0_std": 0.06357885150281098,
    "ZP_5.0_nstars": 30,
    "ZP_6.0": 25.779508960596722,
    "ZP_6.0_std": 0.06472718057714481,
    "ZP_6.0_nstars": 30,
    "ZP_7.0": 25.851914108149213,
    "ZP_7.0_std": 0.06594640100455139,
    "ZP_7.0_nstars": 30,
    "ZP_8.0": 25.89904772872925,
    "ZP_8.0_std": 0.06650479236590319,
    "ZP_8.0_nstars": 30,
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

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )


if __name__ == "__main__":
    print("Calculating latest ZP dictionary")

    new_res, new_errorstack = pipeline.reduce_images(
        dataset=Dataset(ImageBatch()), catch_all_errors=False
    )
