"""
Tests for WINTER reduction
"""
import logging
import unittest

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.97239112854004,
    "ZP_2.0_std": 0.04485376924276352,
    "ZP_2.0_nstars": 20,
    "ZP_3.0": 24.469276428222656,
    "ZP_3.0_std": 0.06575325131416321,
    "ZP_3.0_nstars": 24,
    "ZP_4.0": 24.69712257385254,
    "ZP_4.0_std": 0.05887956544756889,
    "ZP_4.0_nstars": 21,
    "ZP_5.0": 24.742101669311523,
    "ZP_5.0_std": 0.06479335576295853,
    "ZP_5.0_nstars": 22,
    "ZP_6.0": 24.743587493896484,
    "ZP_6.0_std": 0.06024789810180664,
    "ZP_6.0_nstars": 21,
    "ZP_7.0": 24.753494262695312,
    "ZP_7.0_std": 0.06898591667413712,
    "ZP_7.0_nstars": 24,
    "ZP_8.0": 24.76802635192871,
    "ZP_8.0_std": 0.06132316589355469,
    "ZP_8.0_nstars": 21,
    "ZP_AUTO": 24.73688507080078,
    "ZP_AUTO_std": 0.04644254967570305,
    "ZP_AUTO_nstars": 18,
    "SCORMEAN": -0.12958823890710153,
    "SCORMED": -0.1301039138898356,
    "SCORSTD": 1.3122340990145394,
}
expected_dataframe_values = {
    "magpsf": [
        15.375861884283529,
        16.374333222331714,
        15.608207538893376,
        15.067430526438294,
        15.653983533357021,
        13.419238512947432,
        17.011711896680684,
        12.170644129676203,
        13.357850739732427,
    ],
    "magap": [
        15.631678216134478,
        13.205179088230622,
        14.944221301597864,
        13.728325953741994,
        13.330006266001618,
        16.01533818407105,
        13.071709898390827,
        11.974951516927666,
        14.643142758189727,
    ],
}
pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230726"
)

logging.basicConfig(level=logging.DEBUG)


@unittest.skip(
    "WFAU is down, but this test should pass as the values have been"
    "updated in this PR"
)
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

        self.assertEqual(len(candidates_table), 9)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
