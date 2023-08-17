"""
Tests forced photometry with a WIRC image
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.io import open_fits
from mirar.paths import LATEST_SAVE_KEY, ZP_KEY, ZP_STD_KEY, get_output_path
from mirar.pipelines.wirc.blocks import (
    candidate_photometry,
    export_candidates_from_header,
    image_photometry,
)
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_pipeline import WircPipeline
from mirar.processors.utils.header_annotate import HeaderAnnotator, HeaderEditor
from mirar.processors.utils.image_loader import ImageLoader
from mirar.processors.utils.image_saver import ImageSaver
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

NIGHT_NAME = "20210330"

EXPECTED_HEADER_VALUES = {
    "MAGAP": 17.01042620967162,
    "MAGPSF": 16.957491654080247,
    "MAGERRAP": 0.05580721809850005,
    "SIGMAPSF": 0.04793450626938275,
}

EXPECTED_DATAFRAME_VALUES = {
    "magpsf": [16.95749166765038],
    "magap": [17.014130826817812],
    "sigmapsf": [0.04793450628321317],
    "sigmagap": [0.055862250242806104],
    "objectid": ["ZTF21aagppzg"],
}

test_fp_configuration = (
    [
        ImageLoader(
            input_img_dir=test_data_dir.as_posix(),
            input_sub_dir="stack",
            load_image=load_raw_wirc_image,
        ),
        HeaderEditor(
            edit_keys=["TARGRA", "TARGDEC", "TARGNAME"],
            values=[160.643041603707, 34.4374610722322, "ZTF21aagppzg"],
        ),
        HeaderAnnotator(input_keys=["TMC_ZP"], output_key=ZP_KEY),
        HeaderAnnotator(input_keys=["TMC_ZPSD"], output_key=ZP_STD_KEY),
    ]
    + image_photometry
    + [ImageSaver(output_dir_name="photom")]
    + export_candidates_from_header
    + candidate_photometry
)

pipeline = WircPipeline(night=NIGHT_NAME, selected_configurations="test_fp")
pipeline.add_configuration(
    configuration_name="test_fp", configuration=test_fp_configuration
)
pipeline.configure_processors(test_fp_configuration, sub_dir=f"wirc/{NIGHT_NAME}")


class TestForcedPhot(BaseTestCase):
    """
    Class to test WIRC image subtraction pipeline
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Test the pipeline
        """
        self.logger.info("\n\n Testing forced photometry \n\n")

        res, _ = pipeline.reduce_images(
            dataset=Dataset(ImageBatch()), catch_all_errors=False
        )
        self.assertEqual(len(res), 1)

        candidates_table = res[0][0].get_data()
        diff_imgpath = get_output_path(
            base_name=candidates_table.iloc[0][LATEST_SAVE_KEY],
            dir_root="photom",
            sub_dir=NIGHT_NAME,
        )

        # _, header = open_fits(diff_imgpath)
        # for key, value in EXPECTED_HEADER_VALUES.items():
        #     if isinstance(value, float):
        #         self.assertAlmostEqual(value, header[key], places=2)
        #     elif isinstance(value, int):
        #         self.assertEqual(value, header[key])
        #     else:
        #         raise TypeError(
        #             f"Type for value ({type(value)} is neither float not int."
        #         )

        self.assertEqual(len(candidates_table), 1)
        for key, value in EXPECTED_DATAFRAME_VALUES.items():
            if isinstance(value, list):
                for ind, val in enumerate(value):
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )


if __name__ == "__main__":
    print("Calculating latest forced phot dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images(
        dataset=Dataset(ImageBatch()), catch_all_errors=False
    )
    new_candidates_table = new_res[0][0].get_data()
    new_diff_imgpath = get_output_path(
        base_name=new_candidates_table.iloc[0][LATEST_SAVE_KEY],
        dir_root="photom",
        sub_dir="20210330",
    )
    _, new_header = open_fits(new_diff_imgpath)

    NEW_EXP_HEADER = "expected_header_values = { \n"
    for key in EXPECTED_HEADER_VALUES:
        NEW_EXP_HEADER += f'    "{key}": {new_header[key]}, \n'
    NEW_EXP_HEADER += "}"
    print(NEW_EXP_HEADER)

    NEW_EXP_DATAFRAME = "expected_dataframe_values = { \n"
    for key in EXPECTED_DATAFRAME_VALUES:
        NEW_EXP_DATAFRAME += f'    "{key}": {list(new_candidates_table[key])}, \n'
    NEW_EXP_DATAFRAME += "}"
    print(NEW_EXP_DATAFRAME)
