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
    "ZP_4.0": 25.944757945287797,
    "ZP_4.0_std": 0.49857585767841284,
    "ZP_4.0_nstars": 42,
    "ZP_6.0": 26.623309585117156,
    "ZP_6.0_std": 0.4724771413715337,
    "ZP_6.0_nstars": 42,
    "ZP_8.0": 26.992235837754748,
    "ZP_8.0_std": 0.4544249717014124,
    "ZP_8.0_nstars": 42,
    "ZP_10.0": 27.20203013202122,
    "ZP_10.0_std": 0.4537419044965269,
    "ZP_10.0_nstars": 42,
    "ZP_AUTO": 27.03152120235988,
    "ZP_AUTO_std": 0.5298220244029958,
    "ZP_AUTO_nstars": 42,
}

expected_ap_phot = {
    "fluxap2": 67.2952646726674,
    "fluxuncap2": 116.78870886643895,
    "fluxap3": 140.59308642554566,
    "fluxuncap3": 176.46045017437314,
    "magap3": 21.661611289716124,
    "sigmagap3": 1.4624052574158826,
    "fluxap4": 280.95853878914204,
    "fluxuncap4": 234.71042636381088,
    "magap4": 20.909915613483555,
    "sigmagap4": 1.0506128633112535,
    "fluxap5": 284.07456543030924,
    "fluxuncap5": 292.0564971080343,
    "magap5": 20.897940324914178,
    "sigmagap5": 1.2358462219129922,
    "fluxap10": 1096.4122358435654,
    "fluxuncap10": 584.1465261477068,
    "magap10": 19.43158651842342,
    "sigmagap10": 0.7845306178009698,
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
