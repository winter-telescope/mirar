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
    "ZP_2.0": 24.52299118041992,
    "ZP_2.0_std": 0.0678536668419838,
    "ZP_2.0_nstars": 617,
    "ZP_4.0": 25.24053764343262,
    "ZP_4.0_std": 0.04157085344195366,
    "ZP_4.0_nstars": 616,
    "ZP_5.0": 25.32686233520508,
    "ZP_5.0_std": 0.03873402997851372,
    "ZP_5.0_nstars": 615,
    "ZP_8.0": 25.4289379119873,
    "ZP_8.0_std": 0.03831901773810387,
    "ZP_8.0_nstars": 616,
    "ZP_10.0": 25.45285606384277,
    "ZP_10.0_std": 0.03835428133606911,
    "ZP_10.0_nstars": 615,
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
    print("{")
    for key, value in expected_header.items():
        print(f""""{key}": {header[key]},""")
    print("}")
