import unittest
import logging
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.pipelines.summer.summer_files import summer_mask_path, summer_weight_path, \
    sextractor_astrometry_config, sextractor_photometry_config, scamp_path, swarp_path
from winterdrp.processors.utils import ImageSaver
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors import MaskPixels, BiasCalibrator, FlatCalibrator
from winterdrp.processors.csvlog import CSVLog
from winterdrp.paths import core_fields, base_name_key
from winterdrp.pipelines.summer.summer_pipeline import load_raw_summer_image, summer_pixel_scale, \
    summer_astrometric_catalog_generator, summer_photometric_catalog_generator, SummerPipeline
from winterdrp.downloader.get_test_data import get_test_data_dir

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

test_configuration = [
    ImageLoader(
        input_img_dir=test_data_dir,
        input_sub_dir="raw",
        load_image=load_raw_summer_image
    ),
    CSVLog(
        export_keys=[
                        "UTC", 'FIELDID', "FILTERID", "EXPTIME", "OBSTYPE", "RA", "DEC", "TARGTYPE",
                        base_name_key
                    ] + core_fields
    ),
    # DatabaseExporter(
    #     db_name=pipeline_name,
    #     db_table="exposures",
    #     schema_path=get_summer_schema_path("exposures")
    # ),
    MaskPixels(mask_path=summer_mask_path),
    # SplitImage(
    #     buffer_pixels=0,
    #     n_x=1,
    #     n_y=2
    # ),
    # ImageSaver(output_dir_name="rawimages"),
    # DatabaseExporter(
    #     db_name=pipeline_name,
    #     db_table="raw",
    #     schema_path=get_summer_schema_path("raw")
    # ),
    BiasCalibrator(),
    ImageBatcher(split_key="filter"),
    FlatCalibrator(),
    ImageSelector(("OBSTYPE", "SCIENCE")),
    AutoAstrometry(pa=0, inv=True, pixel_scale=summer_pixel_scale),
    Sextractor(
        output_sub_dir="test",
        weight_image=summer_weight_path,
        checkimage_name=None,
        checkimage_type=None,
        **sextractor_astrometry_config
    ),
    Scamp(
        ref_catalog_generator=summer_astrometric_catalog_generator,
        scamp_config_path=scamp_path,
    ),
    Swarp(swarp_config_path=swarp_path, imgpixsize=2400),
    Sextractor(output_sub_dir="test",
               checkimage_type='BACKGROUND_RMS',
               **sextractor_photometry_config),
    PhotCalibrator(ref_catalog_generator=summer_photometric_catalog_generator),
    ImageSaver(output_dir_name="test")
]

expected_zp = {
    "ZP_2.0": 24.381802801640834,
    "ZP_2.0_std": 0.08222868531222044,
    "ZP_2.0_nstars": 30,
    "ZP_3.0": 25.06422489916484,
    "ZP_3.0_std": 0.06870406816807186,
    "ZP_3.0_nstars": 30,
    "ZP_4.0": 25.451230355453497,
    "ZP_4.0_std": 0.06158930864995796,
    "ZP_4.0_nstars": 30,
    "ZP_5.0": 25.659666622034713,
    "ZP_5.0_std": 0.06283258570749543,
    "ZP_5.0_nstars": 30,
    "ZP_6.0": 25.78091181551616,
    "ZP_6.0_std": 0.06342292306060344,
    "ZP_6.0_nstars": 30,
    "ZP_7.0": 25.853334717305508,
    "ZP_7.0_std": 0.06322232852162107,
    "ZP_7.0_nstars": 30,
    "ZP_8.0": 25.90134099756877,
    "ZP_8.0_std": 0.06359669923996372,
    "ZP_8.0_nstars": 30
}

test_config_name = "test"

pipeline = SummerPipeline(night="20220402", selected_configurations=[test_config_name])
pipeline.add_configuration(configuration_name=test_config_name, configuration=test_configuration)


class TestSummerPipeline(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, errorstack = pipeline.reduce_images([[[], []]])
        self.assertEqual(len(errorstack.reports), 0)
        self.assertEqual(len(res[0][0]), 1)

        header = res[0][1][0]

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(f"Type for value ({type(value)} is neither float not int.")


if __name__ == "__main__":

    print("Calculating latest ZP dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images([[[], []]], catch_all_errors=False)

    new_header = new_res[0][1][0]

    new_exp = "expected_zp = { \n"
    for header_key in new_header.keys():
        if "ZP_" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)


