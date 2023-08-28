"""
Tests for WINTER reduction
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.11035919189453,
    "ZP_2.0_std": 0.03795776143670082,
    "ZP_2.0_nstars": 23,
    "ZP_3.0": 24.618452072143555,
    "ZP_3.0_std": 0.02989865280687809,
    "ZP_3.0_nstars": 22,
    "ZP_4.0": 24.83037567138672,
    "ZP_4.0_std": 0.024005483835935593,
    "ZP_4.0_nstars": 18,
    "ZP_5.0": 24.888044357299805,
    "ZP_5.0_std": 0.03897009789943695,
    "ZP_5.0_nstars": 22,
    "ZP_6.0": 24.914031982421875,
    "ZP_6.0_std": 0.04439741000533104,
    "ZP_6.0_nstars": 24,
    "ZP_7.0": 24.910863876342773,
    "ZP_7.0_std": 0.056104876101017,
    "ZP_7.0_nstars": 25,
    "ZP_8.0": 24.896026611328125,
    "ZP_8.0_std": 0.052887771278619766,
    "ZP_8.0_nstars": 24,
    "ZP_AUTO": 24.91974449157715,
    "ZP_AUTO_std": 0.044717829674482346,
    "ZP_AUTO_nstars": 21,
    "SCORMEAN": 0.017224874268454884,
    "SCORMED": 0.014726609007903662,
    "SCORSTD": 1.2294695354674363,
}

expected_dataframe_values = {
    "magpsf": [
        13.45081731465378,
        13.58264248927471,
        16.442835390011275,
        13.651407374516433,
        13.4514507223608,
        13.451014251336805,
    ],
    "magap": [
        11.433320883682018,
        12.039888596367417,
        12.662222810044911,
        11.50348045890572,
        11.518537933186485,
        11.436219364257235,
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

        self.assertEqual(len(candidates_table), 6)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
