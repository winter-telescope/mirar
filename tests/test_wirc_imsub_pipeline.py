"""
Tests for image subtraction with WIRC
"""
import logging
import os

from astropy.io import fits

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.io import open_fits
from mirar.paths import get_output_path
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
    "SCORMEAN": -0.04392849749015427,
    "SCORMED": -0.02331383704510348,
    "SCORSTD": 1.25595271525748,
}

EXPECTED_DATAFRAME_VALUES = {
    "magpsf": [
        19.532174878440024,
        19.37926219106463,
        19.592829160635986,
        17.551198868994298,
        17.197011228688517,
    ],
    "magap": [
        20.391746490573613,
        18.8437585775724,
        19.25446868459551,
        17.74467203323279,
        17.11032933773533,
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
        if os.environ.get("FRITZ_TOKEN", default="") == "":
            raise RuntimeError(
                "No Fritz token. Set environment variable FRITZ_TOKEN to test."
            )
        if os.environ.get("KOWALSKI_TOKEN", default="") == "":
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

        candidates_table = res[0][0].get_data()
        diff_imgpath = get_output_path(
            base_name=candidates_table.iloc[0]["diffimgname"],
            dir_root="subtract",
            sub_dir=NIGHT_NAME,
        )

        _, header = open_fits(diff_imgpath)
        for key, value in EXPECTED_HEADER_VALUES.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        self.assertEqual(len(candidates_table), 5)
        for key, value in EXPECTED_DATAFRAME_VALUES.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )


if __name__ == "__main__":
    print("Calculating latest scorr metrics dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images(
        dataset=Dataset(ImageBatch()), catch_all_errors=False
    )
    new_candidates_table = new_res[0][0].get_data()
    new_diff_imgpath = get_output_path(
        base_name=new_candidates_table.iloc[0]["diffimname"],
        dir_root="subtract",
        sub_dir="20210330",
    )
    _, new_header = open_fits(new_diff_imgpath)

    NEW_EXP_HEADER = "expected_header_values = { \n"
    for header_key in new_header.keys():
        if "SCOR" in header_key:
            NEW_EXP_HEADER += f'    "{header_key}": {new_header[header_key]}, \n'
    NEW_EXP_HEADER += "}"
    print(NEW_EXP_HEADER)

    NEW_EXP_DATAFRAME = "expected_dataframe_values = { \n"
    for key in EXPECTED_DATAFRAME_VALUES:
        NEW_EXP_DATAFRAME += f'    "{key}": {list(new_candidates_table[key])}, \n'
    NEW_EXP_DATAFRAME += "}"
    print(NEW_EXP_DATAFRAME)
