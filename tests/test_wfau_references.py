"""
Tests for getting and making WFAU reference images
"""
import logging

import pandas as pd

from winterdrp.pipelines import get_pipeline
from winterdrp.pipelines.winter.build_references import (
    run_winter_reference_build_pipeline,
    winter_fields_file,
)
from winterdrp.testing import BaseTestCase

logger = logging.getLogger(__name__)

TEST_WINTER_FIELD_ID = 5842
expected_header = {
    "ZP_2.0": 24.12174606323242,
    "ZP_2.0_std": 0.0787130668759346,
    "ZP_2.0_nstars": 273,
    "ZP_4.0": 25.03158950805664,
    "ZP_4.0_std": 0.04695883393287659,
    "ZP_4.0_nstars": 270,
    "ZP_5.0": 25.18747520446777,
    "ZP_5.0_std": 0.0409061461687088,
    "ZP_5.0_nstars": 272,
    "ZP_8.0": 25.36505317687988,
    "ZP_8.0_std": 0.03634443134069443,
    "ZP_8.0_nstars": 273,
    "ZP_10.0": 25.40524482727051,
    "ZP_10.0_std": 0.03677214309573174,
    "ZP_10.0_nstars": 271,
    "COADDS": 12,
    "RA0_0": 282.782828029719,
    "DEC0_0": 45.57267757395213,
    "RA1_0": 282.7921389237596,
    "DEC1_0": 46.51289983153139,
    "RA0_1": 281.6743963051082,
    "DEC0_1": 45.57256743925012,
    "RA1_1": 281.6646948009752,
    "DEC1_1": 46.51278602248236,
    "FIELDID": 5842,
    "SUBDETID": 0,
}

pipeline = get_pipeline(
    instrument="ir_reference_building",
    selected_configurations=["default"],
    night="references",
)


class TestIRReferencePipeline(BaseTestCase):
    """
    Module for testing IR reference building pipeline
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
        Test IR reference building pipeline
        Returns:

        """
        self.logger.info("\n\n Testing IR reference building pipeline \n\n")
        winter_fields = pd.read_csv(winter_fields_file, delim_whitespace=True)

        selected_winter_field = winter_fields[
            winter_fields["ID"] == TEST_WINTER_FIELD_ID
        ].reset_index(drop=True)
        res, _ = run_winter_reference_build_pipeline(
            winter_fields=selected_winter_field,
            nx=1,
            ny=1,
            full_ra_size_deg=0.5,
            full_dec_size_deg=0.6,
        )

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()
        for key, value in expected_header.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )


if __name__ == "__main__":
    print("Calculating WFAU reference image dictionary")

    winter_fields = pd.read_csv(winter_fields_file, delim_whitespace=True)

    selected_winter_field = winter_fields[winter_fields["ID"] == TEST_WINTER_FIELD_ID]
    new_res, new_errorstack = run_winter_reference_build_pipeline(
        winter_fields=selected_winter_field,
        nx=1,
        ny=1,
        full_ra_size_deg=0.5,
        full_dec_size_deg=0.6,
    )

    header = new_res[0][0].get_header()
    print("{")
    for key, value in expected_header.items():
        print(f""""{key}": {header[key]},""")
    print("}")
