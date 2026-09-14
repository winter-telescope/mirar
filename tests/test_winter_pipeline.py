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
        17.228474555700856,
        17.45785163208918,
        17.555629066095435,
        17.32679783287861,
        17.421408291622374,
        17.59762578981333,
        17.683866215375765,
        17.21348601063599,
        17.503135318848273,
        17.54062327531579,
    ],
    "magap": [
        17.42901154178621,
        17.24300685886131,
        16.9354861095572,
        17.396366024618345,
        17.460566708596204,
        17.04155387816826,
        17.00676178707949,
        16.780059358281363,
        16.99542179767127,
        17.55111910454456,
    ],
    # -999.0 is PS1's own sentinel for "no r-band measurement", not a
    # missing/failed crossmatch - the object still matched (see distpsnr1).
    "srmag1": [
        21.1719,
        21.0142,
        -999.0,
        -999.0,
        -999.0,
        -999.0,
        20.5753,
        -999.0,
        22.0182,
        -999.0,
    ],
    "distpsnr1": [
        29.450583335827737,
        26.557587085073617,
        28.50005010290259,
        18.851799554267494,
        26.7737708958688,
        29.150562436876385,
        29.30929330508708,
        29.8918738647234,
        29.49107070373421,
        27.497209944870622,
    ],
}
expected_dataframe_ids = {
    "psobjectid1": [
        172792108188358841,
        172802109896440327,
        172842109802730716,
        172822109253789586,
        172822108630904406,
        172842108838395428,
        172822108916429851,
        172832111076449768,
        172842110466630898,
        172852109512460485,
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
