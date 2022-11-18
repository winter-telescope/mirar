import unittest
from winterdrp.processors.dark import MasterDarkCalibrator
from winterdrp.processors.flat import MasterFlatCalibrator
from winterdrp.processors.sky import MasterSkyCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path, sextractor_astrometry_config, scamp_fp_path, \
    swarp_sp_path
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
import logging
import os
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline
from winterdrp.pipelines.wirc.generator import wirc_astrometric_catalog_generator, wirc_photometric_catalog_generator
from winterdrp.processors.csvlog import CSVLog
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.data import DataSet, ImageBatch

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_2.0": 25.774316787719727,
    "ZP_2.0_std": 0.32469287514686584,
    "ZP_2.0_nstars": 15,
    "ZP_4.0": 27.227283477783203,
    "ZP_4.0_std": 0.09417041391134262,
    "ZP_4.0_nstars": 13,
    "ZP_5.0": 27.59035301208496,
    "ZP_5.0_std": 0.08887326717376709,
    "ZP_5.0_nstars": 13,
    "ZP_8.0": 28.080671310424805,
    "ZP_8.0_std": 0.21738892793655396,
    "ZP_8.0_nstars": 15,
    "ZP_10.0": 28.32463264465332,
    "ZP_10.0_std": 0.08160017430782318,
    "ZP_10.0_nstars": 13,
    "ZP_AUTO": 28.462337493896484,
    "ZP_AUTO_std": 0.18256542086601257,
    "ZP_AUTO_nstars": 15,
}


def get_cal_path(name: str) -> str:
    return os.path.join(test_data_dir, f"wirc/cals/test_{name}.fits")


test_configuration = [
    ImageLoader(
        input_img_dir=test_data_dir,
        input_sub_dir="raw",
        load_image=load_raw_wirc_image
    ),
    CSVLog(
        export_keys=["OBJECT", "FILTER", "UTSHUT", "EXPTIME", "COADDS", "OBSTYPE", "OBSCLASS"],
    ),
    MaskPixels(mask_path=wirc_mask_path),
    ImageSelector(("exptime", "45.0")),
    MasterDarkCalibrator(get_cal_path("dark")),
    ImageDebatcher(),
    ImageSelector(("obsclass", "science")),
    ImageBatcher(split_key="filter"),
    MasterFlatCalibrator(get_cal_path("flat")),
    MasterSkyCalibrator(get_cal_path("sky")),
    ImageSelector(
        ("object", "ZTF21aagppzg"),
        ("filter", "J")
    ),
    AutoAstrometry(catalog="tmc"),
    Sextractor(
        output_sub_dir="postprocess",
        **sextractor_astrometry_config
    ),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
    ),
    Swarp(swarp_config_path=swarp_sp_path),
    Sextractor(
        output_sub_dir="final_sextractor",
        **sextractor_astrometry_config
    ),
    ImageSaver(output_dir_name="final"),
    PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator)
]

pipeline = WircPipeline(night="20210330", selected_configurations="test")
pipeline.add_configuration(configuration_name="test", configuration=test_configuration)


class TestWircPipeline(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing wirc pipeline \n\n")

        res, errorstack = pipeline.reduce_images(DataSet([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res), 1)

        header = res[0][0].get_header()

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=4)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(f"Type for value ({type(value)} is neither float not int.")

