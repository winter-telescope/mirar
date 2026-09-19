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
    "ZP_2.0": 23.97849170284599,
    "ZP_2.0_std": 0.05826721581504818,
    "ZP_2.0_nstars": 60,
    "ZP_3.0": 24.450847552345795,
    "ZP_3.0_std": 0.050819998428748195,
    "ZP_3.0_nstars": 60,
    "ZP_4.0": 24.643778328182464,
    "ZP_4.0_std": 0.04755603507046444,
    "ZP_4.0_nstars": 60,
    "ZP_5.0": 24.729127733225777,
    "ZP_5.0_std": 0.04711641560751583,
    "ZP_5.0_nstars": 60,
    "ZP_6.0": 24.771803826980506,
    "ZP_6.0_std": 0.0467397353333381,
    "ZP_6.0_nstars": 60,
    "ZP_7.0": 24.7951094448713,
    "ZP_7.0_std": 0.04556357052011898,
    "ZP_7.0_nstars": 59,
    "ZP_8.0": 24.83036267022543,
    "ZP_8.0_std": 0.04895778665964604,
    "ZP_8.0_nstars": 59,
    "ZP_AUTO": 24.837867474098875,
    "ZP_AUTO_std": 0.05223561116074162,
    "ZP_AUTO_nstars": 60,
    "ZP_PSF": 24.68237371708096,
    "ZP_PSF_std": 0.05574147047803554,
    "ZP_PSF_nstars": 58,
    "SCORMEAN": -0.13852209047822428,
    "SCORMED": -0.11971798962209092,
    "SCORSTD": 1.366477864242046,
}
expected_dataframe_values = {
    "magpsf": [
        17.459018156318503,
        17.400117283544837,
        17.5704374949843,
        17.421561105370916,
        17.612486746944878,
        17.406527652610777,
        17.53291768741557,
        17.486532379542602,
        17.337293245812347,
        17.59915622384011,
    ],
    "magap": [
        16.980993619145604,
        17.6743745507552,
        17.757116505088504,
        17.404027704001695,
        18.416953527585594,
        17.07614540428232,
        16.98800692894683,
        17.318235944567558,
        17.59918908723602,
        17.128514561978506,
    ],
    "srmag1": [
        None,
        None,
        None,
        18.589899,
        None,
        None,
        None,
        None,
        None,
        21.147499,
    ],
    "distpsnr1": [
        5.303014064214237,
        5.812916238470925,
        4.063396201925486,
        6.538928839304017,
        3.145552403079902,
        3.8404528104622058,
        2.76218077967704,
        5.512099315658212,
        5.759239550127722,
        2.8294753698137916,
    ],
}
expected_dataframe_ids = {
    "psobjectid1": [
        172822108786193609,
        172812107993991984,
        172822110475151072,
        172822110500264488,
        172822108064804690,
        172832109158061138,
        172832109083398985,
        172842109527143353,
        172842108665434558,
        172842108856902396,
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

    def test_pipeline(self):  # pylint: disable=too-many-branches
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

        print("New Results WINTER:")
        print("expected_zp = {")
        for key in expected_zp:
            print(f'    "{key}": {source_table[key]},')
        print("}")

        new_candidates_table = source_table.get_data()

        print("expected_dataframe_values = {")
        for key in expected_dataframe_values:
            print(f'    "{key}": {list(new_candidates_table[key][:10])},')
        print("}")

        print("expected_dataframe_ids = {")
        for key in expected_dataframe_ids:
            print(f'    "{key}": {list(new_candidates_table[key][:10])},')
        print("}")

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

        self.assertEqual(len(candidates_table), 128)
        for key, value in expected_dataframe_values.items():
            for ind, val in enumerate(value):
                if val is None:
                    self.assertIsNone(candidates_table.iloc[ind][key])
                else:
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )

        for key, value in expected_dataframe_ids.items():
            for ind, val in enumerate(value):
                self.assertEqual(candidates_table.iloc[ind][key], val)
