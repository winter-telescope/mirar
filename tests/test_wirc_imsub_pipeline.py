"""
Tests for image subtraction with WIRC
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from winterdrp.data import Dataset, Image, ImageBatch
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.io import open_fits
from winterdrp.paths import get_output_path
from winterdrp.pipelines.wirc.blocks import candidates, image_photometry, subtract
from winterdrp.pipelines.wirc.generator import (
    wirc_reference_image_resampler,
    wirc_reference_psfex,
    wirc_reference_sextractor,
)
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline
from winterdrp.processors.reference import Reference
from winterdrp.processors.utils import ImageLoader, ImageSaver
from winterdrp.processors.utils.header_annotate import HeaderEditor
from winterdrp.references.wirc import WIRCRef
from winterdrp.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

ref_img_directory = test_data_dir.joinpath("wirc/ref")

NIGHT_NAME = "20210330"
expected_candidate_path = Path(test_data_dir).joinpath("wirc/wirc_sources.csv")

expected_res = pd.read_csv(expected_candidate_path)


def test_reference_image_generator(
    image: Image,
    images_directory: str = ref_img_directory,
) -> WIRCRef:
    """
    Generate a test reference image using WIRC

    :param image: science image
    :param images_directory: directiry of reference image
    :return: reference image object
    """
    object_name = image["OBJECT"]
    filter_name = image["FILTER"]
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory,
    )


EXPECTED_HEADER_VALUES = {
    "SCORSTD": 1.081806800432295,
    "SCORMED": -8.757084251543588e-05,
    "SCORMEAN": -0.031172912552408068,
    "MAGAP": 17.104291,
    "MAGPSF": 17.197002,
}

EXPECTED_DATAFRAME_VALUES = {
    "magpsf": [19.319820, 19.242908, 17.197002, 17.565868],
    "magap": [19.302467, 19.122576, 17.104291, 17.917712],
}

test_imsub_configuration = (
    [
        ImageLoader(
            input_img_dir=test_data_dir,
            input_sub_dir="stack",
            load_image=load_raw_wirc_image,
        ),
        Reference(
            ref_image_generator=test_reference_image_generator,
            swarp_resampler=wirc_reference_image_resampler,
            sextractor=wirc_reference_sextractor,
            ref_psfex=wirc_reference_psfex,
        ),
    ]
    + subtract
    + [
        HeaderEditor(
            edit_keys=["TARGRA", "TARGDEC"], values=[160.643041603707, 34.4374610722322]
        )
    ]
    + image_photometry
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
    Class for testing image subtraction with WIRC
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Test the full image subtraction pipeline

        :return: None
        """
        self.logger.info("\n\n Testing wirc imsub pipeline \n\n")

        res, _ = pipeline.reduce_images(
            dataset=Dataset(ImageBatch()), catch_all_errors=False
        )

        self.assertEqual(len(res), 1)

        table = res[0][0].get_data()
        diff_imgpath = get_output_path(
            base_name=table.iloc[0]["diffimname"],
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

        # To update test data, uncomment:
        # table.to_csv(expected_candidate_path, index=False)

        self.assertEqual(len(table), len(expected_res))

        for colname, column in expected_res.items():
            if isinstance(column.iloc[0], (int, np.integer)):
                pd.testing.assert_series_equal(
                    column, table.loc[:, colname], check_dtype=False
                )
            elif isinstance(column.iloc[0], float):
                expected = column.to_numpy(dtype=float)
                res = table.loc[:, colname].to_numpy(dtype=float)
                np.testing.assert_array_almost_equal(expected, res, decimal=2)
            elif colname in ["cutoutScience", "cutoutTemplate", "cutoutDifference"]:
                # Skip cutouts
                pass
            elif isinstance(column.iloc[0], str):
                pd.testing.assert_series_equal(
                    column.astype(str), table.loc[:, colname].astype(str)
                )
            elif "name" in colname:
                # Though path will vary, check base name of images
                self.assertEqual(
                    Path(column.iloc[0]).name, Path(table.loc[0, colname]).name
                )
            else:
                raise ValueError(f"Unexpected column {column}")
