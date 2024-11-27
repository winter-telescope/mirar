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
    "ZP_2.0": 23.993245591021413,
    "ZP_2.0_std": 0.050858550005846714,
    "ZP_2.0_nstars": 63,
    "ZP_3.0": 24.497036449611407,
    "ZP_3.0_std": 0.04728816152282778,
    "ZP_3.0_nstars": 62,
    "ZP_4.0": 24.69026345419674,
    "ZP_4.0_std": 0.045239472510677475,
    "ZP_4.0_nstars": 63,
    "ZP_5.0": 24.783645949516316,
    "ZP_5.0_std": 0.04106679223161631,
    "ZP_5.0_nstars": 63,
    "ZP_6.0": 24.824080937709997,
    "ZP_6.0_std": 0.039621699329724064,
    "ZP_6.0_nstars": 63,
    "ZP_7.0": 24.806677813533444,
    "ZP_7.0_std": 0.04719670888331326,
    "ZP_7.0_nstars": 63,
    "ZP_8.0": 24.821819816703744,
    "ZP_8.0_std": 0.04557932419379016,
    "ZP_8.0_nstars": 63,
    "ZP_AUTO": 24.842581896703418,
    "ZP_AUTO_std": 0.04345809004955868,
    "ZP_AUTO_nstars": 63,
    "ZP_PSF": 24.688038852189116,
    "ZP_PSF_std": 0.046292458825900454,
    "ZP_PSF_nstars": 63,
    "SCORMEAN": -0.11031815756641533,
    "SCORMED": -0.10029606778373226,
    "SCORSTD": 1.237854039383729,
}
expected_dataframe_values = {
    "magpsf": [
        17.249855846497038,
        15.257768513830397,
        17.57158105418562,
        17.626057318913837,
        17.627910755254195,
        17.18914340660828,
        17.517050692176067,
        16.889748567518428,
        17.672502329995737,
        17.100929697191713,
    ],
    "magap": [
        17.96760533096264,
        16.13619912871888,
        16.862310325596575,
        16.98225093924355,
        17.471394905528726,
        17.296224061788823,
        17.00814373242739,
        17.19662557395056,
        17.199305387333204,
        16.930487215194365,
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

        self.assertEqual(len(candidates_table), 153)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
