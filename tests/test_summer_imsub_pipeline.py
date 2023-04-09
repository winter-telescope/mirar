"""
Module to test summer image subtraction pipeline
"""
import logging

from winterdrp.data import Dataset, ImageBatch
from winterdrp.pipelines.summer.summer_pipeline import SummerPipeline
from winterdrp.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_values = {
    "SCORSTD": 1.0490024745298843,
    "SCORMED": 0.0003960102347584023,
    "SCORMEAN": -0.002868667351765475,
}

pipeline = SummerPipeline(night="20220815", selected_configurations=["test_imsub"])


class TestSummerImsubPipeline(BaseTestCase):
    """
    Class to test summer imsub pipeline
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Function to test summer image sub pipeline
        Returns:

        """
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res), 1)

        header = res[0][0].get_header()

        for key, value in expected_values.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )
