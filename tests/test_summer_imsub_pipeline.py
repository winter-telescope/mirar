import unittest
import logging
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector
from winterdrp.paths import base_name_key
from winterdrp.pipelines.summer.summer_pipeline import SummerPipeline
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.pipelines.summer.blocks import subtract
from winterdrp.pipelines.summer.load_summer_image import load_proc_summer_image

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

test_configuration = [
    ImageLoader(
        input_img_dir=test_data_dir,
        input_sub_dir='processed',
        load_image=load_proc_summer_image
    ),
    ImageSelector((base_name_key, "SUMMER_20220816_042349_Camera0.resamp.fits")),
] + subtract


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


