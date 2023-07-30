"""
Module to test SEDMv2 pipeline with "default_transient" configuration
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.pipelines.sedmv2.blocks import process_transient
from mirar.pipelines.sedmv2.load_sedmv2_image import load_sedmv2_mef_image
from mirar.pipelines.sedmv2.sedmv2_pipeline import SEDMv2Pipeline
from mirar.processors.utils import MEFLoader  # MultiExtParser
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_4.0": 25.36942881745319,
    "ZP_4.0_std": 0.14864595223970622,
    "ZP_4.0_nstars": 49,
    "ZP_6.0": 26.162161270312875,
    "ZP_6.0_std": 0.13163613207360314,
    "ZP_6.0_nstars": 49,
    "ZP_8.0": 26.670253060648395,
    "ZP_8.0_std": 0.10943878945493246,
    "ZP_8.0_nstars": 49,
    "ZP_10.0": 27.00897391752905,
    "ZP_10.0_std": 0.08841983371476303,
    "ZP_10.0_nstars": 49,
    "ZP_AUTO": 27.475927448018393,
    "ZP_AUTO_std": 0.1390826192603922,
    "ZP_AUTO_nstars": 42,
}

test_configuration = [
    MEFLoader(
        input_img_dir=test_data_dir,
        input_sub_dir="raw/",
        load_image=load_sedmv2_mef_image,
    ),
] + process_transient

pipeline = SEDMv2Pipeline(night="20230526", selected_configurations="test_transient")
pipeline.add_configuration(
    configuration_name="test_transient", configuration=test_configuration
)


class TestSEDMv2TransientPipeline(BaseTestCase):
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
        function for testing SEDMv2 transient pipeline
        Returns:

        """
        self.logger.info("\n\n Testing SEDMv2 transient pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res), 1)

        header = res[0][0].get_header()

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )


if __name__ == "__main__":
    print("Calculating latest ZP dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images(
        Dataset(ImageBatch()), catch_all_errors=False
    )

    new_header = new_res[0][0].get_header()

    new_exp = "expected_zp = { \n"  # pylint: disable=C0103
    for header_key in new_header.keys():
        if "ZP_" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)
