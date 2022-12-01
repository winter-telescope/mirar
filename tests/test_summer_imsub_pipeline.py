import logging
import unittest

from winterdrp.data import Dataset, ImageBatch
from winterdrp.pipelines.summer.summer_pipeline import SummerPipeline

logger = logging.getLogger(__name__)

expected_values = {
    "SCORSTD": 1.120988782614284,
    "SCORMED": 0.0010565268947477073,
    "SCORMEAN": -0.0027870992375066423,
}

pipeline = SummerPipeline(night="20220815", selected_configurations=["test_imsub"])


class TestSummerPipeline(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, errorstack = pipeline.reduce_images(
            Dataset([ImageBatch()]), catch_all_errors=False
        )

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
