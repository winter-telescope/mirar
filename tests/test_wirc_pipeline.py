"""
Module to test WIRC pipeline
"""
import logging
import os

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.pipelines.wirc.blocks import log, masking, reduction
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_pipeline import WircPipeline
from mirar.processors.dark import MasterDarkCalibrator
from mirar.processors.utils.image_loader import ImageLoader
from mirar.processors.utils.image_selector import ImageSelector
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_6.0": 27.760366439819336,
    "ZP_6.0_std": 0.2747618854045868,
    "ZP_6.0_nstars": 13,
    "ZP_10.0": 28.341157913208008,
    "ZP_10.0_std": 0.09141725301742554,
    "ZP_10.0_nstars": 11,
    "ZP_14.0": 28.417837142944336,
    "ZP_14.0_std": 0.13541148602962494,
    "ZP_14.0_nstars": 13,
    "ZP_18.0": 28.475311279296875,
    "ZP_18.0_std": 0.10811392217874527,
    "ZP_18.0_nstars": 13,
    "ZP_AUTO": 28.48914909362793,
    "ZP_AUTO_std": 0.09673704952001572,
    "ZP_AUTO_nstars": 12,
}


def get_cal_path(name: str) -> str:
    """
    Function to get cal path
    Args:
        name:

    Returns:

    """
    return os.path.join(test_data_dir, f"wirc/cals/test_{name}.fits")


test_configuration = (
    [
        ImageLoader(
            input_img_dir=test_data_dir,
            input_sub_dir="raw",
            load_image=load_raw_wirc_image,
        ),
    ]
    + log
    + masking
    + [ImageSelector(("exptime", "45.0")), MasterDarkCalibrator(get_cal_path("dark"))]
    + reduction
)

pipeline = WircPipeline(night="20210330", selected_configurations="test")
pipeline.add_configuration(configuration_name="test", configuration=test_configuration)


class TestWircPipeline(BaseTestCase):
    """
    Class to Test WIRC Pipeline
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
        function for testing pipeline
        Returns:

        """
        self.logger.info("\n\n Testing wirc pipeline \n\n")

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

    new_exp = "expected_zp = { \n"
    for header_key in new_header.keys():
        if "ZP_" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)
