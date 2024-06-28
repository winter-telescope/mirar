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
    "ZP_2.0": 24.039165434205593,
    "ZP_2.0_std": 0.040682124370149605,
    "ZP_2.0_nstars": 111,
    "ZP_3.0": 24.448507155050052,
    "ZP_3.0_std": 0.0471318743148075,
    "ZP_3.0_nstars": 111,
    "ZP_4.0": 24.635763963930522,
    "ZP_4.0_std": 0.04354799201975238,
    "ZP_4.0_nstars": 111,
    "ZP_5.0": 24.716636555577335,
    "ZP_5.0_std": 0.04169164522178617,
    "ZP_5.0_nstars": 111,
    "ZP_6.0": 24.7431289901993,
    "ZP_6.0_std": 0.041704603416498676,
    "ZP_6.0_nstars": 111,
    "ZP_7.0": 24.76461130331333,
    "ZP_7.0_std": 0.04157117729511644,
    "ZP_7.0_nstars": 111,
    "ZP_8.0": 24.777898334765396,
    "ZP_8.0_std": 0.041735663104911024,
    "ZP_8.0_nstars": 111,
    "ZP_AUTO": 24.783530950220236,
    "ZP_AUTO_std": 0.04301252694189325,
    "ZP_AUTO_nstars": 111,
    "ZP_PSF": 24.662168539775998,
    "ZP_PSF_std": 0.04163357501847084,
    "ZP_PSF_nstars": 110,
    "SCORMEAN": -0.07430693836056637,
    "SCORMED": -0.07059966355274964,
    "SCORSTD": 1.2762548812602645,
}
expected_dataframe_values = {
    "magpsf": [
        16.539472228516754,
        14.816334995321805,
        15.486920360136258,
        17.178854595239663,
        16.357246610808602,
        16.813796280840393,
        17.17955668916707,
        17.208511947122744,
        14.795105119473721,
        16.686258217055684,
    ],
    "magap": [
        16.320732009563713,
        15.413317709017791,
        16.223866004796673,
        18.095085591114476,
        16.959889773064127,
        16.88041688053796,
        19.116076621089814,
        16.984351826663126,
        15.198017344346496,
        16.06515983909224,
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

        self.assertEqual(len(candidates_table), 1274)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )
