import unittest
import logging
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.pipelines.summer.summer_files import summer_mask_path, summer_weight_path, \
    sextractor_astrometry_config, sextractor_photometry_config, scamp_path, swarp_config_path
from winterdrp.processors.utils import ImageSaver
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors import MaskPixels, BiasCalibrator, FlatCalibrator
from winterdrp.processors.csvlog import CSVLog
from winterdrp.paths import core_fields, base_name_key
from winterdrp.pipelines.summer.summer_pipeline import load_proc_summer_image, summer_pixel_scale, \
    summer_reference_sextractor, summer_reference_image_generator, summer_reference_psfex, summer_reference_image_resampler, SummerPipeline
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.processors.reference import Reference
from winterdrp.processors.astromatic import PSFex
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare, default_summer_catalog_purifier

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

test_configuration = [
            ImageLoader(
                input_img_dir=test_data_dir,
                input_sub_dir='processed',
                load_image=load_proc_summer_image
            ),
            ImageBatcher(split_key=base_name_key),
            ImageSelector(('OBSTYPE', 'SCIENCE')),
            # ImageSelector(('FILTER', ['u'])),
            ImageSelector((base_name_key, "SUMMER_20220816_042349_Camera0.resamp.fits")),
            Reference(ref_image_generator=summer_reference_image_generator,
                      ref_psfex=summer_reference_psfex,
                      ref_sextractor=summer_reference_sextractor,
                      ref_swarp_resampler=summer_reference_image_resampler),
            Sextractor(config_path='winterdrp/pipelines/summer/summer_imsub_files/config/photomCat.sex',
                       parameter_path='winterdrp/pipelines/summer/summer_imsub_files/config/photom.param',
                       filter_path='winterdrp/pipelines/summer/summer_imsub_files/config/default.conv',
                       starnnw_path='winterdrp/pipelines/summer/summer_imsub_files/config/default.nnw',
                       output_sub_dir='subtract',
                       cache=False,
                       write_regions_file=True),
            PSFex(config_path='winterdrp/pipelines/summer/summer_imsub_files/config/photom.psfex',
                  output_sub_dir="subtract",
                  norm_fits=True),
            ImageSaver(output_dir_name='ref'),
            ZOGYPrepare(output_sub_dir="subtract", sci_zp_header_key='ZP_AUTO', catalog_purifier=default_summer_catalog_purifier),
            ZOGY(output_sub_dir="subtract"),
            # DetectCandidates(output_sub_dir="subtract",
            #                  cand_det_sextractor_config='winterdrp/pipelines/summer/summer_imsub_files/config/photomCat.sex',
            #                  cand_det_sextractor_nnw='winterdrp/pipelines/summer/summer_imsub_files/config/default.nnw',
            #                  cand_det_sextractor_filter='winterdrp/pipelines/summer/summer_imsub_files/config/default.conv',
            #                  cand_det_sextractor_params='winterdrp/pipelines/summer/summer_imsub_files/config/Scorr.param'),
            # RegionsWriter(output_dir_name='candidates'),
            # PSFPhotometry(),
            # AperturePhotometry(aper_diameters=[8, 40], cutout_size_aper_phot=100, bkg_in_diameters=[25, 90],
            #                    bkg_out_diameters=[40, 100], col_suffix_list=['', 'big']),
            # DataframeWriter(output_dir_name='candidates'),
        ]

expected_values = {
    'SCORSTD': 1.120988782614284,
    'SCORMED': 0.0010565268947477073,
    'SCORMEAN': -0.0027870992375066423
}

test_config_name = "test_imsub"

pipeline = SummerPipeline(night="20220815", selected_configurations=[test_config_name])
pipeline.add_configuration(configuration_name=test_config_name, configuration=test_configuration)


class TestSummerPipeline(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, errorstack = pipeline.reduce_images([[[], []]], catch_all_errors=False)

        self.assertEqual(len(res[0][0]), 1)

        header = res[0][1][0]

        for key, value in expected_values.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(f"Type for value ({type(value)} is neither float not int.")


if __name__ == "__main__":

    print("Calculating latest scorr metrics dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images([[[], []]], catch_all_errors=False)

    new_header = new_res[0][1][0]

    new_exp = "expected_values = { \n"
    for header_key in new_header.keys():
        if "SCOR" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)


