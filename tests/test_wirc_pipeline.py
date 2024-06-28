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
    "ZP_6.0": 28.756793975830078,
    "ZP_6.0_std": 0.07175969332456589,
    "ZP_6.0_nstars": 10,
    "ZP_10.0": 29.221546173095703,
    "ZP_10.0_std": 0.08415297418832779,
    "ZP_10.0_nstars": 10,
    "ZP_14.0": 29.307411193847656,
    "ZP_14.0_std": 0.12999354302883148,
    "ZP_14.0_nstars": 11,
    "ZP_18.0": 29.351654052734375,
    "ZP_18.0_std": 0.10415466874837875,
    "ZP_18.0_nstars": 11,
    "ZP_AUTO": 29.370691299438477,
    "ZP_AUTO_std": 0.0895678773522377,
    "ZP_AUTO_nstars": 11,
    "DETMAG95": 20.52475357055664,
    "DETMAG50": 19.24155044555664,
    "DETMAG05": 16.64298439025879,
    "MAGLIM": 21.068127497605154,
}

expected_table = {
    "ra": 160.64316666666662,
    "dec": 34.437416666666664,
    "xpos": 689.7351634112279,
    "ypos": 1658.0823528760743,
    "fluxap": 35875.90509640997,
    "fluxuncap": 687.4760147074132,
    "magap": 17.983684133207536,
    "sigmagap": 0.09195371412681115,
    "fluxapbig": 74376.7707068217,
    "fluxuncapbig": 3028.8886771963785,
    "magapbig": 17.192098003889477,
    "sigmagapbig": 0.09989157687293666,
    "psfflux": 39665.9926674688,
    "psffluxunc": 417.21051569081976,
    "chipsf": 0.9663571880988399,
    "xshift": 1,
    "yshift": 0,
    "magpsf": 17.874645481282734,
    "sigmapsf": 0.09029330822879632,
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

        # Cleanup
        output_dir = get_output_dir("wirc/20210330")
        if output_dir.exists():
            shutil.rmtree(output_dir)

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
                self.assertAlmostEqual(value, header[key], places=2)
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
            new_exp += f'    "{header_key}": {row[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_table.items():
            if isinstance(value, float):
                ratio = value / row[key]
                self.assertAlmostEqual(ratio, 1, delta=0.005)
            elif isinstance(value, int):
                self.assertEqual(value, row[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)}) is neither float not int."
                )
