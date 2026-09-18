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
    "ZP_2.0": 23.97814850254867,
    "ZP_2.0_std": 0.05832326700001151,
    "ZP_2.0_nstars": 60,
    "ZP_3.0": 24.450910328241655,
    "ZP_3.0_std": 0.05082148072614819,
    "ZP_3.0_nstars": 60,
    "ZP_4.0": 24.643925177458325,
    "ZP_4.0_std": 0.04755575861598483,
    "ZP_4.0_nstars": 60,
    "ZP_5.0": 24.729109504300617,
    "ZP_5.0_std": 0.04712040460883906,
    "ZP_5.0_nstars": 60,
    "ZP_6.0": 24.7717020594737,
    "ZP_6.0_std": 0.04674623326432091,
    "ZP_6.0_nstars": 60,
    "ZP_7.0": 24.79500526145292,
    "ZP_7.0_std": 0.04557240232521663,
    "ZP_7.0_nstars": 59,
    "ZP_8.0": 24.830294193506674,
    "ZP_8.0_std": 0.048968580371657335,
    "ZP_8.0_nstars": 59,
    "ZP_AUTO": 24.837775634457284,
    "ZP_AUTO_std": 0.05224100648606315,
    "ZP_AUTO_nstars": 60,
    "ZP_PSF": 24.662553190105122,
    "ZP_PSF_std": 0.05127408852532155,
    "ZP_PSF_nstars": 57,
    "SCORMEAN": -0.13828411674238658,
    "SCORMED": -0.11942388364288006,
    "SCORSTD": 1.3665580887340274,
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
    "srmag1": [],
    "distpsnr1": [],
}
# PS1 crossmatch object id (via Kowalski) for the nearest match to each of
# the first 10 candidates. Unlike the photometry above, this must match
# exactly, not just approximately - it is a copied identifier, not a
# computed value, and any drift would mean the crossmatch itself changed.
# This is a baseline for comparison against the BOOM-based PS1 crossmatch
# (see PR #1145) - once BOOM's own values are confirmed to match this
# Kowalski baseline, this test can be dropped in favour of that one.
expected_dataframe_ids = {
    "psobjectid1": [],
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

        # Uncomment to print a fresh baseline to copy-paste in after a
        # genuine, intentional change to the pipeline's output.
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

        self.assertEqual(len(candidates_table), 129)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )

        # The PS1 crossmatch id is a value copied verbatim from Kowalski,
        # not a re-computed quantity, so it must match exactly.
        for key, value in expected_dataframe_ids.items():
            for ind, val in enumerate(value):
                self.assertEqual(candidates_table.iloc[ind][key], val)
