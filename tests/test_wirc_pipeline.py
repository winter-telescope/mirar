"""
Module to test WIRC pipeline
"""

import logging
import shutil
from pathlib import Path

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import get_output_dir
from mirar.pipelines.wirc.blocks import log, masking, test
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_pipeline import WircPipeline
from mirar.processors.dark import MasterDarkCalibrator
from mirar.processors.utils.image_loader import ImageLoader
from mirar.processors.utils.image_selector import ImageSelector
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_header = {
    "ZP_6.0": 28.654869079589844,
    "ZP_6.0_std": 0.27347737550735474,
    "ZP_6.0_nstars": 13,
    "ZP_10.0": 29.23480987548828,
    "ZP_10.0_std": 0.09129349142313004,
    "ZP_10.0_nstars": 11,
    "ZP_14.0": 29.311601638793945,
    "ZP_14.0_std": 0.13507914543151855,
    "ZP_14.0_nstars": 13,
    "ZP_18.0": 29.368980407714844,
    "ZP_18.0_std": 0.10783938318490982,
    "ZP_18.0_nstars": 13,
    "ZP_AUTO": 29.382888793945312,
    "ZP_AUTO_std": 0.0966394692659378,
    "ZP_AUTO_nstars": 12,
    "DETMAG95": 20.556614112854000,
    "DETMAG50": 19.223386764526400,
    "DETMAG05": 16.460859966278100,
    "MAGLIM": 21.0444307896674,
}

expected_table = {
    "ra": 160.64316666666662,
    "dec": 34.437416666666664,
    "xpos": 689.7351634112279,
    "ypos": 1658.0823528760743,
    "fluxap": 36445.87475183021,
    "fluxuncap": 700.6073201070337,
    "magap": 17.978775054092406,
    "sigmagap": 0.09919814358921,
    "fluxapbig": 70467.29198086978,
    "fluxuncapbig": 3087.953360586301,
    "magapbig": 17.26292704469811,
    "sigmagapbig": 0.108024192790888,
    "psfflux": 40020.58541113894,
    "psffluxunc": 426.734305980583,
    "chipsf": 0.9351921013740141,
    "xshift": 1,
    "yshift": 0,
    "magpsf": 17.87718740725503,
    "sigmapsf": 0.09766545011532148,
}


def get_test_dark_path(_) -> Path:
    """
    Function to get cal path
    Args:

    Returns:
        :return: Path to calibration file
    """
    return Path(test_data_dir).joinpath("wirc/cals/test_dark.fits")


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
    + [
        ImageSelector(("exptime", "45.0")),
        MasterDarkCalibrator(master_image_path_generator=get_test_dark_path),
    ]
    + test
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

        # Cleanup
        output_dir = get_output_dir("wirc/20210330")
        shutil.rmtree(output_dir)

        header = res[0][0].get_metadata()

        # Check the values from the header, which are derived from the stack

        print("New Results:")
        new_exp = "expected_header = { \n"
        for header_key in expected_header:
            new_exp += f'    "{header_key}": {header[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_header.items():
            if isinstance(value, float):
                ratio = value / header[key]
                self.assertAlmostEqual(ratio, 1, delta=0.005)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)}) is neither float not int."
                )

        src_table = res[0][0].get_data()

        self.assertEqual(len(src_table), 1)

        # Check the values from the forced photometry, applied to the difference image

        row = src_table.iloc[0]

        print("New Results:")
        new_exp = "expected_table = { \n"
        for header_key in expected_table:
            new_exp += f'    "{header_key.lower()}": {row[header_key.lower()]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_table.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, row[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, row[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)}) is neither float not int."
                )
