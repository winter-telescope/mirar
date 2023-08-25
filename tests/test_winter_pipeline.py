"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.344806671142578,
    "ZP_2.0_std": 0.09197966009378433,
    "ZP_2.0_nstars": 34,
    "ZP_3.0": 24.81342887878418,
    "ZP_3.0_std": 0.12821219861507416,
    "ZP_3.0_nstars": 36,
    "ZP_4.0": 24.981557846069336,
    "ZP_4.0_std": 0.1145731657743454,
    "ZP_4.0_nstars": 36,
    "ZP_5.0": 25.02898406982422,
    "ZP_5.0_std": 0.11246983706951141,
    "ZP_5.0_nstars": 36,
    "ZP_6.0": 25.04960060119629,
    "ZP_6.0_std": 0.11517446488142014,
    "ZP_6.0_nstars": 36,
    "ZP_7.0": 25.06537437438965,
    "ZP_7.0_std": 0.11556794494390488,
    "ZP_7.0_nstars": 36,
    "ZP_8.0": 25.075380325317383,
    "ZP_8.0_std": 0.1130020022392273,
    "ZP_8.0_nstars": 36,
    "ZP_AUTO": 25.061182022094727,
    "ZP_AUTO_std": 0.11833330243825912,
    "ZP_AUTO_nstars": 36,
    "SCORMEAN": -0.006769005587283054,
    "SCORMED": -0.016401447652707915,
    "SCORSTD": 1.1897194960119415,
}

expected_dataframe_values = {
    "magpsf": [
        17.45071330284148,
        17.30534842251346,
        17.310026927117633,
        16.82961340224697,
        17.127554097741523,
        16.80151788735043,
        16.793798532434806,
        17.159071102789405,
        16.626942595496722,
        15.885022332662135,
        16.517233281991267,
        17.151807714647745,
    ],
    "magap": [
        13.213044754544207,
        14.910611357578492,
        15.312879230359306,
        14.871036613037449,
        14.02161902260337,
        13.832161416781771,
        14.381031692225502,
        14.909421958453404,
        17.85001937426495,
        10.983572804661817,
        11.212310335581147,
        12.625275407985267,
    ],
}

pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230726"
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

        # Expect one dataset, for one different sub-boards
        self.assertEqual(len(res[0]), 1)

        source_table = res[0][0]

        # # Uncomment to print new expected ZP dict
        print("New Results WINTER:")
        new_exp = "expected_zp = { \n"
        for header_key in source_table.get_metadata():
            if header_key in expected_zp:
                new_exp += f'    "{header_key}": {source_table[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        new_candidates_table = source_table.get_data()

        new_exp_dataframe = "expected_dataframe_values = { \n"
        for key in expected_dataframe_values:
            new_exp_dataframe += f'    "{key}": {list(new_candidates_table[key])}, \n'
        new_exp_dataframe += "}"

        print(new_exp_dataframe)

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, source_table[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, source_table[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        candidates_table = source_table.get_data()

        self.assertEqual(len(candidates_table), 40)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
