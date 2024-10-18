"""
Tests for WINTER reduction
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.paths import get_output_dir
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 23.98942717700761,
    "ZP_2.0_std": 0.04344578250755102,
    "ZP_2.0_nstars": 112,
    "ZP_3.0": 24.477112607639857,
    "ZP_3.0_std": 0.042484225250743156,
    "ZP_3.0_nstars": 111,
    "ZP_4.0": 24.654207750819698,
    "ZP_4.0_std": 0.04439331556411262,
    "ZP_4.0_nstars": 112,
    "ZP_5.0": 24.7221475232846,
    "ZP_5.0_std": 0.04519185139825257,
    "ZP_5.0_nstars": 112,
    "ZP_6.0": 24.756044837456937,
    "ZP_6.0_std": 0.044349808138425866,
    "ZP_6.0_nstars": 112,
    "ZP_7.0": 24.772752692282218,
    "ZP_7.0_std": 0.04482551006672038,
    "ZP_7.0_nstars": 112,
    "ZP_8.0": 24.770982276168485,
    "ZP_8.0_std": 0.04579597308653655,
    "ZP_8.0_nstars": 112,
    "ZP_AUTO": 24.77890676740926,
    "ZP_AUTO_std": 0.04744931807399898,
    "ZP_AUTO_nstars": 112,
    "ZP_PSF": 24.642527455478035,
    "ZP_PSF_std": 0.0399356259304595,
    "ZP_PSF_nstars": 111,
    "SCORMEAN": -0.07475470995979283,
    "SCORMED": -0.0707062139750039,
    "SCORSTD": 1.2678094266462752,
}
expected_dataframe_values = {
    "magpsf": [
        17.15442302684727,
        15.185488779330399,
        16.98839381259969,
        17.142272198245134,
        16.087920413156862,
        16.414188801000417,
        16.22731390503941,
        15.803972171134065,
        17.226034715962143,
        15.215588861635071,
    ],
    "magap": [
        17.450959465495885,
        15.913252980944264,
        16.53822890545202,
        16.896747559635436,
        16.450315549584033,
        17.519977695006027,
        16.649724469850867,
        16.00580147048224,
        17.813509037544115,
        15.518791251926956,
    ],
}

pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230726"
)

logging.basicConfig(level=logging.DEBUG)


# @unittest.skip(
#     "WFAU is down"
# )
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

        # Cleanup - delete ouptut dir
        output_dir = get_output_dir(dir_root="winter/20230726")
        shutil.rmtree(output_dir)

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
            new_exp_dataframe += (
                f'    "{key}": {list(new_candidates_table[key][:10])}, \n'
            )
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

        self.assertEqual(len(candidates_table), 61)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
