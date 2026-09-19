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
        17.458812542360995,
        17.399940714760763,
        17.570775919978885,
        17.421680982116335,
        17.612249713190142,
        17.406640865693245,
        17.532911856690518,
        17.48649965211238,
        17.337396644716204,
        17.59887440605349,
    ],
    "magap": [
        16.98031087845887,
        17.680104004152902,
        17.758904889497096,
        17.42562922682663,
        18.410150871089456,
        17.07635726369595,
        16.987159652173418,
        17.308115455430396,
        17.597206293933148,
        17.124283337139047,
    ],
    "srmag1": [
        None,
        None,
        21.9872,
        None,
        21.0019,
        21.820299,
        None,
        21.697001,
        None,
        20.887898999999997,
    ],
    "distpsnr1": [
        None,
        5.926316031879033,
        14.273739577173421,
        7.365720721493199,
        15.563648925556727,
        15.980351600904857,
        9.966652427453882,
        5.761397381177185,
        10.996305734150413,
        14.308060315020628,
    ],
}
expected_dataframe_ids = {
    "psobjectid1": [
        None,
        1.728121079937006e17,
        1.7281211047136765e17,
        1.728221104649579e17,
        1.7282210814785398e17,
        1.7282210918875494e17,
        1.7284210913936003e17,
        1.728421095578917e17,
        1.7283210868259917e17,
        1.7283210888185677e17,
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
