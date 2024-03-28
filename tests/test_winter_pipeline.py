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
    "ZP_2.0": 23.959992146021143,
    "ZP_2.0_std": 0.04468593502379257,
    "ZP_2.0_nstars": 118,
    "ZP_3.0": 24.488647137271013,
    "ZP_3.0_std": 0.043678000162985946,
    "ZP_3.0_nstars": 118,
    "ZP_4.0": 24.67403631422544,
    "ZP_4.0_std": 0.040829710238539714,
    "ZP_4.0_nstars": 118,
    "ZP_5.0": 24.736553404147642,
    "ZP_5.0_std": 0.040354725833759444,
    "ZP_5.0_nstars": 118,
    "ZP_6.0": 24.77188170392595,
    "ZP_6.0_std": 0.039765732833753016,
    "ZP_6.0_nstars": 118,
    "ZP_7.0": 24.782397765144424,
    "ZP_7.0_std": 0.04069653155959119,
    "ZP_7.0_nstars": 118,
    "ZP_8.0": 24.799413106226112,
    "ZP_8.0_std": 0.040876147446689284,
    "ZP_8.0_nstars": 118,
    "ZP_AUTO": 24.787539557156656,
    "ZP_AUTO_std": 0.04316430194985392,
    "ZP_AUTO_nstars": 118,
    "ZP_PSF": 24.680094141593663,
    "ZP_PSF_std": 0.03875121055985624,
    "ZP_PSF_nstars": 117,
    "SCORMEAN": -0.07678255815161257,
    "SCORMED": -0.07307016879473982,
    "SCORSTD": 1.2836771224604269,
}
expected_dataframe_values = {
    "magpsf": [
        16.39444674612076,
        15.405810813687662,
        17.082173218426902,
        17.082294915192634,
        17.403596735252293,
        17.160933991469214,
        14.66480196304238,
        14.657302634953442,
        16.816141837382137,
        16.25453113498711,
    ],
    "magap": [
        16.48067858723591,
        15.915934341171843,
        18.53190378870282,
        18.264227223646593,
        16.264918566468864,
        17.122984755361113,
        15.45802381805369,
        15.202372311912047,
        16.160387918403462,
        16.456876102361825,
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

        self.assertEqual(len(candidates_table), 1295)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
