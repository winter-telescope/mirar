"""
Tests forced photometry with a WIRC image
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import ZP_KEY, ZP_STD_KEY, get_output_dir
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
    "TARGNAME": "ZTF21aagppzg",
}

EXPECTED_DATAFRAME_VALUES = {
    "magpsf": [16.95749166765038],
    "magap": [17.014130826817812],
    "sigmapsf": [0.04793450628321317],
    "sigmagap": [0.055862250242806104],
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

        res, _, _ = pipeline.reduce_images(
            dataset=Dataset(ImageBatch()), catch_all_errors=False
        )
        self.assertEqual(len(res), 1)
        # Cleanup
        output_dir = get_output_dir(f"wirc/{NIGHT_NAME}")
        shutil.rmtree(output_dir)

        new = res[0][0]

        source_table = new.get_data()
        metadata = new.get_metadata()

        print("Got the following values:")

        new_header = "expected_header_values = { \n"
        for key in EXPECTED_HEADER_VALUES:
            new_header += f'    "{key}": {metadata[key]}, \n'
        new_header += "}"
        print(new_header)

        new_src_table = "expected_dataframe_values = { \n"
        for key in EXPECTED_DATAFRAME_VALUES:
            new_src_table += f'    "{key}": {list(source_table[key])}, \n'
        new_src_table += "}"
        print(new_src_table)

        for key, value in EXPECTED_HEADER_VALUES.items():
            if isinstance(value, str):
                self.assertEqual(value, metadata[key])
            else:
                self.assertAlmostEqual(value, metadata[key], delta=0.05)

        self.assertEqual(len(source_table), 1)

        for key, value in EXPECTED_DATAFRAME_VALUES.items():
            for i, row in source_table.iterrows():
                self.assertAlmostEqual(row[key], value[i], delta=0.05)
