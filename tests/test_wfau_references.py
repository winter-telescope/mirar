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
    "ZP_2.0": 24.511962890625,
    "ZP_2.0_std": 0.06453035771846771,
    "ZP_2.0_nstars": 260,
    "ZP_4.0": 25.23262596130371,
    "ZP_4.0_std": 0.04047837853431702,
    "ZP_4.0_nstars": 259,
    "ZP_5.0": 25.32016944885254,
    "ZP_5.0_std": 0.0390668548643589,
    "ZP_5.0_nstars": 259,
    "ZP_8.0": 25.42309188842773,
    "ZP_8.0_std": 0.0383964516222477,
    "ZP_8.0_nstars": 257,
    "ZP_10.0": 25.44612503051758,
    "ZP_10.0_std": 0.03797683864831924,
    "ZP_10.0_nstars": 256,
    "COADDS": 4,
    "RA0_0": 282.0557287433173,
    "DEC0_0": 45.35327585919234,
    "RA1_0": 282.0582425176955,
    "DEC1_0": 45.82890665679109,
    "RA0_1": 281.3771483296388,
    "DEC0_1": 45.35301704115277,
    "RA1_1": 281.3738911418988,
    "DEC1_1": 45.82864350802002,
    "FIELDID": 5842,
    "SUBDETID": 0,
}

pipeline = get_pipeline(
    instrument="winter",
    selected_configurations=["refbuild"],
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
        self.logger.info("\n\n Testing WINTER reference building pipeline \n\n")
        winter_fields = pd.read_csv(winter_fields_file, delim_whitespace=True)

        selected_winter_field = winter_fields[
            winter_fields["ID"] == TEST_WINTER_FIELD_ID
        ].reset_index(drop=True)
        res, _ = run_winter_reference_build_pipeline(
            winter_fields=selected_winter_field, only_this_subdet_id=0
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
        winter_fields=selected_winter_field, only_this_subdet_id=0
    )

    header = new_res[0][0].get_header()
    print("{")
    for key, value in expected_header.items():
        print(f""""{key}": {header[key]},""")
    print("}")
