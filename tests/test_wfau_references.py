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
    "ZP": 25.412062499999998,
    "ZPSTD": 0.003714472329619354,
    "COADDS": 4,
    "RA0_0": 282.0717301991488,
    "DEC0_0": 45.341829671758134,
    "RA1_0": 282.0746173634903,
    "DEC1_0": 45.81745868212202,
    "RA0_1": 281.3931283796939,
    "DEC0_1": 45.341830141404145,
    "RA1_1": 281.3902425647145,
    "DEC1_1": 45.81745915963218,
    "FIELDID": 5842,
    "SUBDETID": 0,
}

pipeline = get_pipeline(
    instrument="winter",
    selected_configurations=["refbuild"],
    night="references",
)


# @unittest.skip("WFAU is down")
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
            subdet_id=0, field_id=TEST_WINTER_FIELD_ID, catch_all_errors=False
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
