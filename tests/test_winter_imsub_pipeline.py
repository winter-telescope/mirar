"""
Tests for WINTER image subtraction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.io import open_raw_image
from mirar.paths import DIFF_IMG_KEY, get_output_path
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

NIGHT_NAME = "20230726"

EXPECTED_HEADER_VALUES = {
    "SCORMEAN": -0.13625905938349483,
    "SCORMED": -0.13512726465303718,
    "SCORSTD": 1.3021071256299737,
}
EXPECTED_DATAFRAME_VALUES = {
    "magpsf": [12.173105476828795],
    "magap": [12.062389584233465],
}


pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test_imsub"], night="20230726"
)

logging.basicConfig(level=logging.DEBUG)


class TestWinterImsubPipeline(BaseTestCase):
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
        self.logger.info("\n\n Testing winter imsub pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Expect two datasets, for two different sub-boards
        self.assertEqual(len(res), 1)

        candidates_table = res[0][0].get_data()
        diff_imgpath = get_output_path(
            base_name=candidates_table.iloc[0][DIFF_IMG_KEY],
            dir_root="diffs",
            sub_dir=NIGHT_NAME,
        )

        diff_image = open_raw_image(diff_imgpath)

        print("New Results WINTER IMSUB:")
        new_exp = "EXPECTED_HEADER_VALUES = { \n"
        for header_key in diff_image.header.keys():
            if header_key in EXPECTED_HEADER_VALUES:
                new_exp += f'    "{header_key}": {diff_image[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        new_df = "EXPECTED_DATAFRAME_VALUES = { \n"
        for key in candidates_table.keys():
            if key in EXPECTED_DATAFRAME_VALUES:
                new_df += f'    "{key}": {candidates_table[key].tolist()}, \n'
        new_df += "}"
        print(new_df)

        for key, value in EXPECTED_HEADER_VALUES.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, diff_image[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, diff_image[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        self.assertEqual(len(candidates_table), 1)
        for key, value in EXPECTED_DATAFRAME_VALUES.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
