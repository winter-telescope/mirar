import unittest
from winterdrp.processors.dark import DarkCalibrator
from winterdrp.processors.flat import SkyFlatCalibrator
from winterdrp.processors.sky import NightSkyMedianCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
from winterdrp.paths import coadd_key, proc_history_key
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.processors.dark import DarkCalibrator
from winterdrp.processors.flat import SkyFlatCalibrator
from winterdrp.processors.sky import NightSkyMedianCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path, sextractor_astrometry_config, scamp_fp_path, \
    swarp_sp_path
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.catalog import Gaia2Mass
from winterdrp.downloader.caltech import download_via_ssh
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
from winterdrp.paths import coadd_key, proc_history_key
from winterdrp.processors.csvlog import CSVLog
import logging
import os
from winterdrp.pipelines.wirc.wirc_pipeline import load_raw_wirc_image, WircPipeline, \
    wirc_astrometric_catalog_generator, wirc_photometric_catalog_generator
from winterdrp.processors.csvlog import CSVLog


logger = logging.getLogger(__name__)

test_data_dir = "/Users/robertstein/Code/winterdrp_staterpack/"

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

test_pipeline = [
    ImageLoader(
        input_sub_dir="raw",
        load_image=load_raw_wirc_image
    ),
    CSVLog(
        export_keys=["OBJECT", "FILTER", "UTSHUT", "EXPTIME", "COADDS", "OBSTYPE", "OBSCLASS"]
    ),
    MaskPixels(mask_path=wirc_mask_path),
    ImageSelector(("exptime", "45.0")),
    DarkCalibrator(lambda x: os.path.join(test_data_dir, "test_dark.fits")),
    ImageDebatcher(),
    ImageSelector(("obsclass", "science")),
    ImageBatcher(split_key="filter"),
    SkyFlatCalibrator(lambda x: os.path.join(test_data_dir, "test_flat.fits")),
    NightSkyMedianCalibrator(lambda x: os.path.join(test_data_dir, "test_sky.fits")),
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

pipeline = WircPipeline(pipeline_configuration=test_pipeline, night="20210330")


class TestWircPipeline(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing wirc pipeline \n\n")

        res, errorstack = pipeline.reduce_images([[[], []]])

        self.assertEqual(len(res), 1)

        header = res[0][1][0]

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key])
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(f"Type for value ({type(value)} is neither float not int.")


if __name__ == "__main__":

    print("Calculating latest ZP dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images([[[], []]])

    new_header = new_res[0][1][0]

    new_exp = "expected_zp = { \n"
    for header_key in new_header.keys():
        if "ZP_" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)

