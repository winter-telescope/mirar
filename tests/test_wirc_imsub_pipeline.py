"""
Tests for image subtraction with WIRC
"""

import logging
import os
import shutil

from astropy.io import fits

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import get_output_dir
from mirar.pipelines.wirc.blocks import candidates, subtract
from mirar.pipelines.wirc.generator import (
    wirc_reference_image_resampler,
    wirc_reference_psfex,
    wirc_reference_sextractor,
)
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_pipeline import WircPipeline
from mirar.processors.reference import ProcessReference
from mirar.processors.utils.image_loader import ImageLoader
from mirar.processors.utils.image_saver import ImageSaver
from mirar.references import WIRCRef
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

ref_img_directory = test_data_dir.joinpath("wirc/ref")

NIGHT_NAME = "20210330"


def reference_image_test_generator(
    header: fits.header,
    images_directory: str = ref_img_directory,
):
    """
    Function to generate reference image for testing
    Args:
        header: image header
        images_directory: reference image directory

    Returns:
        Reference Image
    """
    object_name = header["OBJECT"]
    filter_name = header["FILTER"]
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory,
    )


EXPECTED_HEADER_VALUES = {
    "SCORMEAN": 0.049450589357321745,
    "SCORMED": 0.06992433914552572,
    "SCORSTD": 1.2498783381487255,
}
EXPECTED_DATAFRAME_VALUES = {
    "magpsf": [
        19.37980722187377,
        19.173776427588024,
        17.16839618364596,
        17.593767904409027,
    ],
    "magap": [
        19.68843701280298,
        19.417567597329914,
        17.312072330263064,
        17.902854210820834,
    ],
}

test_imsub_configuration = (
    [
        ImageLoader(
            input_img_dir=test_data_dir.as_posix(),
            input_sub_dir="stack",
            load_image=load_raw_wirc_image,
        ),
        ProcessReference(
            ref_image_generator=reference_image_test_generator,
            swarp_resampler=wirc_reference_image_resampler,
            sextractor=wirc_reference_sextractor,
            ref_psfex=wirc_reference_psfex,
        ),
    ]
    + subtract
    + [ImageSaver(output_dir_name="subtract")]
    + candidates
)

pipeline = WircPipeline(night=NIGHT_NAME, selected_configurations="test_imsub")
pipeline.add_configuration(
    configuration_name="test_imsub", configuration=test_imsub_configuration
)
pipeline.configure_processors(test_imsub_configuration)


class TestWircImsubPipeline(BaseTestCase):
    """
    Class to test WIRC image subtraction pipeline
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_tokens()

    def check_tokens(self):
        """If required tokens do not exist raise an error.
        :raises RuntimeError: If missing token.
        """
        if os.getenv("FRITZ_TOKEN", default="") == "":
            raise RuntimeError(
                "No Fritz token. Set environment variable FRITZ_TOKEN to test."
            )
        if os.getenv("KOWALSKI_TOKEN", default="") == "":
            raise RuntimeError(
                "No Kowalski token. Set environment variable KOWALSKI_TOKEN to test."
            )

    def test_pipeline(self):
        """
        Test the wirc imsub pipeline
        """
        self.logger.info("\n\n Testing wirc imsub pipeline \n\n")

        res, _ = pipeline.reduce_images(
            dataset=Dataset(ImageBatch()), catch_all_errors=False
        )

        self.assertEqual(len(res), 1)

        source_table = res[0][0]
        candidates_table = source_table.get_data()

        print("New Results WIRC-imsub:")
        new_exp_header = "EXPECTED_HEADER_VALUES = { \n"
        for key, _ in EXPECTED_HEADER_VALUES.items():
            new_exp_header += f'"{key}": {source_table[key]},\n'
        new_exp_header += "}"
        print(new_exp_header)

        new_exp_df = "EXPECTED_DATAFRAME_VALUES = { \n"
        for key, _ in EXPECTED_DATAFRAME_VALUES.items():
            new_exp_df += f'"{key}": {list(candidates_table[key])},\n'
        new_exp_df += "}"
        print(new_exp_df)

        for key, value in EXPECTED_HEADER_VALUES.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, source_table[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, source_table[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        self.assertEqual(len(candidates_table), 4)
        for key, value in EXPECTED_DATAFRAME_VALUES.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )

        # Cleanup
        output_dir = get_output_dir("wirc/20210330")
        shutil.rmtree(output_dir)
