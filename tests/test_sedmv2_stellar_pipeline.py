"""
Module to test SEDMv2 pipeline with "default_stellar" configuration
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import get_output_dir
from mirar.pipelines.sedmv2.blocks import image_photometry, process_stellar
from mirar.pipelines.sedmv2.load_sedmv2_image import load_sedmv2_mef_image
from mirar.pipelines.sedmv2.sedmv2_pipeline import SEDMv2Pipeline
from mirar.processors.utils import MEFLoader
from mirar.testing import (
    BaseTestCase,
    load_expected_values,
    save_expected_values,
    update_expected_values,
)

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

EXPECTED = load_expected_values("sedmv2_stellar")
expected_zp = EXPECTED.get("expected_zp", {})
expected_ap_phot = EXPECTED.get("expected_ap_phot", {})


test_configuration = (
    [
        MEFLoader(
            input_img_dir=test_data_dir,
            input_sub_dir="",
            load_image=load_sedmv2_mef_image,
        ),
    ]
    + process_stellar
    + image_photometry
)

pipeline = SEDMv2Pipeline(night="20230328", selected_configurations="test_stellar")
pipeline.add_configuration(
    configuration_name="test_stellar", configuration=test_configuration
)


class TestSEDMv2StellarPipeline(BaseTestCase):
    """
    Class to test SEDMv2 transient pipeline
    """

    def setUp(self):
        """
        Test set up
        Returns:

        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        function for testing SEDMv2 stellar pipeline
        Returns:

        """
        self.check_expected_values("sedmv2_stellar", EXPECTED)
        self.logger.info("\n\n Testing SEDMv2 stellar pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res), 29)

        new = res[0][0]

        header = new.get_metadata()

        if update_expected_values():
            save_expected_values(
                "sedmv2_stellar",
                {
                    "expected_zp": {x: header[x] for x in expected_zp},
                    "expected_ap_phot": {
                        x: new.get_data().iloc[0][x] for x in expected_ap_phot
                    },
                },
            )
            return

        for key, value in expected_zp.items():  # pylint: disable=R0801
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):  # pylint: disable=R0801
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        source_table = new.get_data()

        assert len(source_table) == 1

        for _, row in source_table.iterrows():
            for key, value in expected_ap_phot.items():
                if isinstance(value, float):
                    self.assertAlmostEqual(value, row[key], places=2)
                elif isinstance(value, int):
                    self.assertEqual(value, row[key])
                else:
                    raise TypeError(
                        f"Type for value ({type(value)} is neither float not int."
                    )

        # Cleanup
        output_dir = get_output_dir("sedmv2/20230328")
        shutil.rmtree(output_dir)


if __name__ == "__main__":  # pylint: disable=R0801
    print("Calculating latest ZP dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images(
        Dataset(ImageBatch()), catch_all_errors=False
    )  # pylint: disable=R0801

    new_header = new_res[0][0].get_header()

    new_exp = "expected_zp = { \n"  # pylint: disable=C0103
    for header_key in new_header.keys():
        if "ZP_" in header_key:  # pylint: disable=R0801
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)
