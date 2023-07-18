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
    "ZP_4.0": 25.37294329566108,
    "ZP_4.0_std": 0.15462107846416773,
    "ZP_4.0_nstars": 45,
    "ZP_6.0": 26.173137675679527,
    "ZP_6.0_std": 0.13210938173513917,
    "ZP_6.0_nstars": 45,
    "ZP_8.0": 26.68036031751845,
    "ZP_8.0_std": 0.10794668404233344,
    "ZP_8.0_nstars": 45,
    "ZP_10.0": 27.019062413397897,
    "ZP_10.0_std": 0.08550068701694684,
    "ZP_10.0_nstars": 45,
    "ZP_AUTO": 27.501423087867526,
    "ZP_AUTO_std": 0.059945737356448876,
    "ZP_AUTO_nstars": 37,
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
