"""
Module to test WIRC pipeline
"""
import logging
import os

from winterdrp.data import Dataset, ImageBatch
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.pipelines.wirc.generator import (
    wirc_astrometric_catalog_generator,
    wirc_photometric_catalog_generator,
)
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.wirc_files import (
    scamp_fp_path,
    sextractor_astrometry_config,
    swarp_sp_path,
    wirc_mask_path,
)
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline
from winterdrp.processors.astromatic import Scamp, Sextractor, Swarp
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.csvlog import CSVLog
from winterdrp.processors.dark import MasterDarkCalibrator
from winterdrp.processors.flat import MasterFlatCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.sky import MasterSkyCalibrator
from winterdrp.processors.utils import ImageSaver
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import (
    ImageBatcher,
    ImageDebatcher,
    ImageSelector,
)
from winterdrp.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_2.0": 25.908477783203125,
    "ZP_2.0_std": 0.11761965602636337,
    "ZP_2.0_nstars": 12,
    "ZP_4.0": 27.235090255737305,
    "ZP_4.0_std": 0.1022496446967125,
    "ZP_4.0_nstars": 12,
    "ZP_5.0": 27.5942440032959,
    "ZP_5.0_std": 0.09162138402462006,
    "ZP_5.0_nstars": 12,
    "ZP_8.0": 28.16297149658203,
    "ZP_8.0_std": 0.07523621618747711,
    "ZP_8.0_nstars": 12,
    "ZP_10.0": 28.324859619140625,
    "ZP_10.0_std": 0.07975760102272034,
    "ZP_10.0_nstars": 12,
    "ZP_AUTO": 28.434249877929688,
    "ZP_AUTO_std": 0.14610908925533295,
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


test_configuration = [
    ImageLoader(
        input_img_dir=test_data_dir, input_sub_dir="raw", load_image=load_raw_wirc_image
    ),
    CSVLog(
        export_keys=[
            "OBJECT",
            "FILTER",
            "UTSHUT",
            "EXPTIME",
            "COADDS",
            "OBSTYPE",
            "OBSCLASS",
        ],
    ),
    MaskPixels(mask_path=wirc_mask_path),
    ImageSelector(("exptime", "45.0")),
    MasterDarkCalibrator(get_cal_path("dark")),
    ImageDebatcher(),
    ImageSelector(("obsclass", "science")),
    ImageBatcher(split_key="filter"),
    MasterFlatCalibrator(get_cal_path("flat")),
    MasterSkyCalibrator(get_cal_path("sky")),
    ImageSelector(("object", "ZTF21aagppzg"), ("filter", "J")),
    AutoAstrometry(catalog="tmc"),
    Sextractor(output_sub_dir="postprocess", **sextractor_astrometry_config),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
    ),
    Swarp(swarp_config_path=swarp_sp_path),
    Sextractor(output_sub_dir="final_sextractor", **sextractor_astrometry_config),
    PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator),
    ImageSaver(output_dir_name="final"),
]

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
