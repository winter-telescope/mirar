import unittest
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
from winterdrp.pipelines.wirc_imsub import WircImsubPipeline
from winterdrp.processors.csvlog import CSVLog
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.paths import base_name_key
from winterdrp.processors.reference import Reference
from winterdrp.processors.astromatic import PSFex
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare
from winterdrp.processors.candidates.candidate_detector import DetectCandidates
from winterdrp.pipelines.wirc_imsub.wirc_imsub_pipeline import wirc_reference_psfex, wirc_reference_image_generator, wirc_reference_image_resampler, wirc_reference_sextractor
from winterdrp.pipelines.wirc.wirc_pipeline import load_raw_wirc_image
from astropy.io import fits
from winterdrp.references import WIRCRef

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

ref_img_directory = os.path.join(test_data_dir, 'wirc/ref')
print(ref_img_directory)

def reference_image_generator(
        header: fits.header,
        images_directory: str = ref_img_directory,
):
    object_name = header['OBJECT']
    filter_name = header['FILTER']
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory
    )


expected_values = {
    'SCORSTD': 1.081806800432295,
    'SCORMED': -8.757084251543588e-05,
    'SCORMEAN': -0.031172912552408068
}

test_imsub_configuration = [
    ImageLoader(
        input_img_dir=test_data_dir,
        input_sub_dir="raw",
        load_image=load_raw_wirc_image
    ),
    # ImageBatcher(split_key='UTSHUT'),
    ImageSelector((base_name_key, "ZTF21aagppzg_J_stack_1_20210330.fits")),
    Reference(
        ref_image_generator=reference_image_generator,
        ref_swarp_resampler=wirc_reference_image_resampler,
        ref_sextractor=wirc_reference_sextractor,
        ref_psfex=wirc_reference_psfex
    ),
    # Swarp(),
    Sextractor(config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photomCat.sex',
               parameter_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photom.param',
               filter_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.conv',
               starnnw_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.nnw',
               output_sub_dir='subtract',
               cache=False),
    PSFex(config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photom.psfex',
          output_sub_dir="subtract",
          norm_fits=True),
    ZOGYPrepare(output_sub_dir="subtract"),
    ZOGY(output_sub_dir="subtract")
]

pipeline = WircImsubPipeline(night="20210330", selected_configurations="test_imsub")
pipeline.add_configuration(configuration_name="test_imsub", configuration=test_imsub_configuration)


class TestWircImsubcPipeline(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing wirc imsub pipeline \n\n")

        res, errorstack = pipeline.reduce_images([[[], []]], catch_all_errors=False)

        self.assertEqual(len(res), 1)

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

