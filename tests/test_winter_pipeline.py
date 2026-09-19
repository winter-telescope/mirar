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
        17.459017622216685,
        17.40011826200793,
        17.570436571358826,
        17.421561942578407,
        17.61248707663176,
        17.406526631458902,
        17.53291798648281,
        17.486532049210364,
        17.337293020234085,
        17.59915505301026,
    ],
    "magap": [
        16.980995907134535,
        17.674375545171017,
        17.757124422457988,
        17.40403009975467,
        18.417016172122636,
        17.076147011665505,
        16.98800702585568,
        17.31823587105105,
        17.599194438964943,
        17.128515024174842,
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
        20.887899,
    ],
    "distpsnr1": [
        None,
        5.902478470149945,
        14.268246014324504,
        7.372606579351607,
        15.507445453655324,
        16.001125164778156,
        9.956455829452688,
        5.777884683967926,
        11.06528550272264,
        14.27671213374158,
    ],
    "sgscore1": [
        None,
        0.7536249756813049,
        0.03354166820645332,
        0.003958333283662796,
        0.6446011662483215,
        0.0062500000931322575,
        0.9768750071525574,
        0.4070476293563843,
        0.06437499821186066,
        0.20056547224521637,
    ],
}
expected_dataframe_ids = {
    "psobjectid1": [
        None,
        172812107993700607,
        172812110471367657,
        172822110464957890,
        172822108147853975,
        172822109188754944,
        172842109139360037,
        172842109557891697,
        172832108682599174,
        172832108881856782,
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
