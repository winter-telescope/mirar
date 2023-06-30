"""
Tests for getting and making WFAU reference images
"""
import logging

from mirar.pipelines import get_pipeline
from mirar.pipelines.winter.build_references import run_winter_reference_build_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

TEST_WINTER_FIELD_ID = 5842
expected_header = {
    "ZP_2.0": 24.51775550842285,
    "ZP_2.0_std": 0.06673479825258255,
    "ZP_2.0_nstars": 522,
    "ZP_4.0": 25.237825393676758,
    "ZP_4.0_std": 0.041513219475746155,
    "ZP_4.0_nstars": 522,
    "ZP_5.0": 25.324222564697266,
    "ZP_5.0_std": 0.03828829526901245,
    "ZP_5.0_nstars": 520,
    "ZP_8.0": 25.427032470703125,
    "ZP_8.0_std": 0.038392938673496246,
    "ZP_8.0_nstars": 521,
    "ZP_10.0": 25.450834274291992,
    "ZP_10.0_std": 0.0386686846613884,
    "ZP_10.0_nstars": 521,
    "COADDS": 4,
    "RA0_0": 282.05572874331733,
    "DEC0_0": 45.353275859192344,
    "RA1_0": 282.0582425176955,
    "DEC1_0": 45.82890665679109,
    "RA0_1": 281.37714832963877,
    "DEC0_1": 45.353017041152775,
    "RA1_1": 281.37389114189875,
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

        res, _ = run_winter_reference_build_pipeline(
            subdet_id=0, field_id=TEST_WINTER_FIELD_ID
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

    new_res, new_errorstack = run_winter_reference_build_pipeline(
        field_id=TEST_WINTER_FIELD_ID, subdet_id=0
    )

    header = new_res[0][0].get_header()
    print("expected_header = {")
    for key, value in expected_header.items():
        print(f""""{key}": {header[key]},""")
    print("}")
