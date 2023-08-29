"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

<<<<<<< HEAD
expected_zp = {
    "ZP_2.0": 23.98996925354004,
    "ZP_2.0_std": 0.03564070165157318,
    "ZP_2.0_nstars": 30,
    "ZP_3.0": 24.477495193481445,
    "ZP_3.0_std": 0.04281936213374138,
    "ZP_3.0_nstars": 31,
    "ZP_4.0": 24.63871955871582,
    "ZP_4.0_std": 0.04711835831403732,
    "ZP_4.0_nstars": 32,
    "ZP_5.0": 24.724422454833984,
    "ZP_5.0_std": 0.0430118553340435,
    "ZP_5.0_nstars": 30,
    "ZP_6.0": 24.78891944885254,
    "ZP_6.0_std": 0.03930766507983208,
    "ZP_6.0_nstars": 27,
    "ZP_7.0": 24.79821014404297,
    "ZP_7.0_std": 0.04727752134203911,
    "ZP_7.0_nstars": 30,
    "ZP_8.0": 24.81728744506836,
    "ZP_8.0_std": 0.05584121495485306,
    "ZP_8.0_nstars": 32,
    "ZP_AUTO": 24.794212341308594,
    "ZP_AUTO_std": 0.043879956007003784,
    "ZP_AUTO_nstars": 29,
    "SCORMEAN": -0.1324857697402311,
    "SCORMED": -0.1305649672434809,
    "SCORSTD": 1.3039875768136526,
}
expected_dataframe_values = {
    "magpsf": [12.110608868454948],
    "magap": [12.020240804368768],
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

        self.assertEqual(len(candidates_table), 1)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
