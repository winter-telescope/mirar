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
    "ZP_2.0": 24.005294438914305,
    "ZP_2.0_std": 0.06005590150517328,
    "ZP_2.0_nstars": 59,
    "ZP_3.0": 24.47122127292256,
    "ZP_3.0_std": 0.05191250056994636,
    "ZP_3.0_nstars": 59,
    "ZP_4.0": 24.660134047602313,
    "ZP_4.0_std": 0.04895237732433349,
    "ZP_4.0_nstars": 59,
    "ZP_5.0": 24.73844478378757,
    "ZP_5.0_std": 0.04843877831790673,
    "ZP_5.0_nstars": 59,
    "ZP_6.0": 24.780598421278032,
    "ZP_6.0_std": 0.04778940744176772,
    "ZP_6.0_nstars": 59,
    "ZP_7.0": 24.79966483902917,
    "ZP_7.0_std": 0.045979386039383334,
    "ZP_7.0_nstars": 57,
    "ZP_8.0": 24.836408870088746,
    "ZP_8.0_std": 0.049043100282806104,
    "ZP_8.0_nstars": 56,
    "ZP_AUTO": 24.84092339714826,
    "ZP_AUTO_std": 0.05281705601550301,
    "ZP_AUTO_nstars": 59,
    "ZP_PSF": 24.677819766367843,
    "ZP_PSF_std": 0.05191139701034343,
    "ZP_PSF_nstars": 57,
    "SCORMEAN": -0.13718873361312706,
    "SCORMED": -0.11876823466321687,
    "SCORSTD": 1.3589516850692291,
}
expected_dataframe_values = {
    "magpsf": [
        17.39802916470696,
        17.5095656886213,
        17.40702084758793,
        17.576428667161935,
        17.390064760561224,
        17.503951609537477,
        17.45941660520178,
        17.47893024399312,
        17.325970378802964,
        17.217715061027427,
    ],
    "magap": [
        17.460308769278605,
        16.866293511023073,
        17.406897608738145,
        18.176968259952798,
        17.13565906248416,
        16.923007184191857,
        17.479606721546585,
        17.2121344211263,
        17.592897139280986,
        16.86732637398696,
    ],
    "srmag1": [None, None, None, 21.0019, 21.820299, None, None, 21.697001, None, None],
    "distpsnr1": [
        5.728989073771944,
        None,
        7.414251184284056,
        15.556980114051607,
        16.13806502009323,
        10.114689758646128,
        0.70704665975367,
        5.844587757261037,
        11.001578870846398,
        29.471923214945615,
    ],
    "sgscore1": [
        0.7536249756813049,
        None,
        0.003958333283662796,
        0.6446011662483215,
        0.0062500000931322575,
        0.9768750071525574,
        0.06363095343112946,
        0.4070476293563843,
        0.06437499821186066,
        0.06363095343112946,
    ],
}
expected_dataframe_ids = {
    "psobjectid1": [
        172812107993700607,
        None,
        172822110464957890,
        172822108147853975,
        172822109188754944,
        172842109139360037,
        172832111005989103,
        172842109557891697,
        172832108682599174,
        172832111005989103,
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

        print(f"expected_n_candidates = {len(new_candidates_table)}")

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
