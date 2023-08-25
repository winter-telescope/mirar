"""
Module to test SEDMv2 pipeline with "default_stellar" configuration
"""
import logging

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.pipelines.sedmv2.blocks import image_photometry, process_stellar
from mirar.pipelines.sedmv2.load_sedmv2_image import load_sedmv2_mef_image
from mirar.pipelines.sedmv2.sedmv2_pipeline import SEDMv2Pipeline
from mirar.processors.utils import MEFLoader  # MultiExtParser
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_4.0": 25.94123243758373,
    "ZP_4.0_std": 0.49819926577461493,
    "ZP_4.0_nstars": 39,
    "ZP_6.0": 26.614818280400986,
    "ZP_6.0_std": 0.4728077007528148,
    "ZP_6.0_nstars": 39,
    "ZP_8.0": 26.98211547740056,
    "ZP_8.0_std": 0.4582100267087655,
    "ZP_8.0_nstars": 39,
    "ZP_10.0": 27.189519198354084,
    "ZP_10.0_std": 0.4580645022400429,
    "ZP_10.0_nstars": 39,
    "ZP_AUTO": 27.009442904897835,
    "ZP_AUTO_std": 0.5465585617804521,
    "ZP_AUTO_nstars": 39,
}


expected_ap_phot = {
    "fluxap2": -10.440357798826348,
    "fluxuncap2": 205.54719177533997,
    "fluxap3": 185.81936010251917,
    "fluxuncap3": 310.73759182222466,
    "magap3": 21.336715504439553,
    "sigmagap3": 1.8965330652977908,
    "fluxap4": 288.53853734268324,
    "fluxuncap4": 412.89523151752974,
    "magap4": 20.85893334004581,
    "sigmagap4": 1.6473637846943183,
    "fluxap5": 164.94609445952742,
    "fluxuncap5": 513.6055760247314,
    "magap5": 21.466087812602073,
    "sigmagap5": 3.425448357780633,
    "fluxap10": 804.3668543646124,
    "fluxuncap10": 1027.5301316831465,
    "magap10": 19.745807489218134,
    "sigmagap10": 1.4910821986784057,
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

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        self.assertEqual(len(res), 29)

        new = res[0][0]

        header = new.get_metadata()

        new_exp = "expected_zp = { \n"  # pylint: disable=C0103
        for header_key in expected_zp.keys():
            new_exp += f'    "{header_key}": {header[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        source_table = new.get_data()

        assert len(source_table) == 1

        new_row = "expected_ap_phot = { \n"  # pylint: disable=C0103
        for header_key in expected_ap_phot.keys():
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


if __name__ == "__main__":
    print("Calculating latest ZP dictionary")

    # Code to generate updated ZP dict of the results change

    new_res, new_errorstack = pipeline.reduce_images(
        Dataset(ImageBatch()), catch_all_errors=False
    )

    new_header = new_res[0][0].get_header()

    new_exp = "expected_zp = { \n"  # pylint: disable=C0103
    for header_key in new_header.keys():
        if "ZP_" in header_key:
            new_exp += f'    "{header_key}": {new_header[header_key]}, \n'
    new_exp += "}"
    print(new_exp)
