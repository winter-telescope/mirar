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
# PS1 crossmatch object id (via BOOM's PS1_DR2 catalog) for the nearest match
# to each of the first 10 candidates. Unlike the photometry above, this must
# match exactly, not just approximately - it is a copied identifier, not a
# computed value, and any drift would mean the crossmatch itself changed.
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


def _print_new_expected(var_name: str, items) -> None:
    """
    Print a freshly-measured dict in the same literal form used to pin
    expected_zp/expected_dataframe_values/expected_dataframe_ids above,
    so a new baseline can be copy-pasted in after a genuine, intentional
    change to the pipeline's output.
    """
    body = "".join(f'    "{key}": {value}, \n' for key, value in items)
    print(f"{var_name} = {{ \n{body}}}")


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
        _print_new_expected(
            "expected_zp",
            (
                (key, source_table[key])
                for key in source_table.get_metadata()
                if key in expected_zp
            ),
        )

        new_candidates_table = source_table.get_data()

        _print_new_expected(
            "expected_dataframe_values",
            (
                (key, list(new_candidates_table[key][:10]))
                for key in expected_dataframe_values
            ),
        )
        _print_new_expected(
            "expected_dataframe_ids",
            (
                (key, list(new_candidates_table[key][:10]))
                for key in expected_dataframe_ids
            ),
        )

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

        self.assertEqual(len(candidates_table), 108)
        for key, value in expected_dataframe_values.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )

        # The PS1 crossmatch id is a value copied verbatim from the BOOM
        # response, not a re-computed quantity, so it must match exactly.
        # This guards against the PS1/PS1SGSc/PS1STRM query merge (or any
        # future BOOM schema change) silently matching a different object.
        for key, value in expected_dataframe_ids.items():
            for ind, val in enumerate(value):
                self.assertEqual(candidates_table.iloc[ind][key], val)
