import unittest
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
import logging
import os
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.paths import base_name_key
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline
from astropy.io import fits
from winterdrp.references import WIRCRef
from winterdrp.processors.reference import Reference
from winterdrp.pipelines.wirc.blocks import subtract
from winterdrp.pipelines.wirc.generator import wirc_reference_image_resampler, wirc_reference_sextractor, \
    wirc_reference_psfex
try:
    from winterdrp.data import Dataset, ImageBatch
except ImportError:
    pass


logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO)

test_data_dir = get_test_data_dir()

ref_img_directory = os.path.join(test_data_dir, 'wirc/ref')


def test_reference_image_generator(
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
        input_sub_dir="stack",
        load_image=load_raw_wirc_image
    ),
    Reference(
        ref_image_generator=test_reference_image_generator,
        ref_swarp_resampler=wirc_reference_image_resampler,
        ref_sextractor=wirc_reference_sextractor,
        ref_psfex=wirc_reference_psfex
    ),
] + subtract

pipeline = WircPipeline(night="20210330", selected_configurations="test_imsub")
pipeline.add_configuration(configuration_name="test_imsub", configuration=test_imsub_configuration)


class TestWircImsubPipeline(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing wirc imsub pipeline \n\n")

        try:
            res, errorstack = pipeline.reduce_images(dataset=Dataset(ImageBatch()), catch_all_errors=False)
        except NameError:
            pipeline.reduce_images([[[], []]], catch_all_errors=False)

        self.assertEqual(len(res), 1)

        header = res[0][0].get_header()

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

    new_res, new_errorstack = pipeline.reduce_images(catch_all_errors=False)

    new_header = new_res[0][1][0]

    new_exp = "expected_values = { \n"
    for header_key in new_header.keys():
        if "SCOR" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)

