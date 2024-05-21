"""
Module to test SEDMv2 pipeline with "default_stellar" configuration
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import get_output_dir
from mirar.pipelines.sedmv2.blocks import image_photometry, process_stellar
from mirar.pipelines.sedmv2.load_sedmv2_image import load_sedmv2_mef_image
from mirar.pipelines.sedmv2.sedmv2_pipeline import SEDMv2Pipeline
from mirar.processors.utils import MEFLoader
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_4.0": 25.914995160146862,
    "ZP_4.0_std": 0.4942314614792097,
    "ZP_4.0_nstars": 39,
    "ZP_6.0": 26.58560238463573,
    "ZP_6.0_std": 0.46460089687204953,
    "ZP_6.0_nstars": 39,
    "ZP_8.0": 26.953794458262127,
    "ZP_8.0_std": 0.4464802951237832,
    "ZP_8.0_nstars": 39,
    "ZP_10.0": 27.16306114531297,
    "ZP_10.0_std": 0.44522379930341205,
    "ZP_10.0_nstars": 39,
    "ZP_AUTO": 26.97540065904275,
    "ZP_AUTO_std": 0.5444036865100269,
    "ZP_AUTO_nstars": 39,
}

expected_ap_phot = {
    "fluxap2": -60.76028299697255,
    "fluxuncap2": 57.95930827143862,
    "fluxap3": 26.5731727246716,
    "fluxuncap3": 89.86650958054098,
    "magap3": 23.414251883010365,
    "sigmagap3": 3.7128065432283983,
    "fluxap4": 124.0908628521507,
    "fluxuncap4": 118.88546080295447,
    "magap4": 21.74101089811672,
    "sigmagap4": 1.174224669453573,
    "fluxap5": 93.9257000535545,
    "fluxuncap5": 146.761956351453,
    "magap5": 22.043399307442723,
    "sigmagap5": 1.782072957351011,
    "fluxap10": 5.162133327971844,
    "fluxuncap10": 297.2813813457539,
    "magap10": 25.193287365736992,
    "sigmagap10": 62.54387222395732,
}

test_configuration = (
    [
        MEFLoader(
            input_img_dir=test_data_dir,
            input_sub_dir="",
            load_image=load_sedmv2_mef_image,
        ),
    ]
    + process_stellar
    + image_photometry
)

pipeline = SEDMv2Pipeline(night="20230328", selected_configurations="test_stellar")
pipeline.add_configuration(
    configuration_name="test_stellar", configuration=test_configuration
)


class TestSEDMv2StellarPipeline(BaseTestCase):
    """
    Class to test SEDMv2 transient pipeline
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
        function for testing SEDMv2 stellar pipeline
        Returns:

        """
        self.logger.info("\n\n Testing SEDMv2 stellar pipeline \n\n")

        res, _, _ = pipeline.reduce_images(
            Dataset([ImageBatch()]), catch_all_errors=False
        )

        self.assertEqual(len(res), 29)

        new = res[0][0]

        header = new.get_metadata()

        new_exp = "expected_zp = { \n"  # pylint: disable=C0103, W0621
        for header_key in expected_zp.keys():  # pylint: disable=C0201, W0621
            new_exp += f'    "{header_key}": {header[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_zp.items():  # pylint: disable=R0801
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):  # pylint: disable=R0801
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        source_table = new.get_data()

        assert len(source_table) == 1

        new_row = "expected_ap_phot = { \n"  # pylint: disable=C0103
        for header_key in expected_ap_phot.keys():  # pylint: disable=C0201
            new_row += f'    "{header_key}": {source_table.iloc[0][header_key]}, \n'
        new_row += "}"
        print(new_row)

        for _, row in source_table.iterrows():
            for key, value in expected_ap_phot.items():
                if isinstance(value, float):
                    self.assertAlmostEqual(value, row[key], places=2)
                elif isinstance(value, int):
                    self.assertEqual(value, row[key])
                else:
                    raise TypeError(
                        f"Type for value ({type(value)} is neither float not int."
                    )

        # Cleanup
        output_dir = get_output_dir("sedmv2/20230328")
        shutil.rmtree(output_dir)


if __name__ == "__main__":  # pylint: disable=R0801
    print("Calculating latest ZP dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images(
        Dataset(ImageBatch()), catch_all_errors=False
    )  # pylint: disable=R0801

    new_header = new_res[0][0].get_header()

    new_exp = "expected_zp = { \n"  # pylint: disable=C0103
    for header_key in new_header.keys():
        if "ZP_" in header_key:  # pylint: disable=R0801
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)
