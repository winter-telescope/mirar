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
from mirar.processors.utils import MEFLoader  # MultiExtParser
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

expected_zp = {
    "ZP_4.0": 25.940887427559876,
    "ZP_4.0_std": 0.49824802025845366,
    "ZP_4.0_nstars": 39,
    "ZP_6.0": 26.614804158685146,
    "ZP_6.0_std": 0.4727982155483261,
    "ZP_6.0_nstars": 39,
    "ZP_8.0": 26.982064199066162,
    "ZP_8.0_std": 0.4582201430387493,
    "ZP_8.0_nstars": 39,
    "ZP_10.0": 27.189493400241165,
    "ZP_10.0_std": 0.45806133571891794,
    "ZP_10.0_nstars": 39,
    "ZP_AUTO": 27.008659106875687,
    "ZP_AUTO_std": 0.5468613179756767,
    "ZP_AUTO_nstars": 39,
}

expected_ap_phot = {
    "fluxap2": -13.850341348155112,
    "fluxuncap2": 205.58084553347814,
    "fluxap3": 161.299057848526,
    "fluxuncap3": 310.77711898765494,
    "magap3": 21.489579530194685,
    "sigmagap3": 2.1626931885857243,
    "fluxap4": 302.3792599800752,
    "fluxuncap4": 412.9893558727741,
    "magap4": 20.807279107272198,
    "sigmagap4": 1.5808578175981052,
    "fluxap5": 103.20045314743703,
    "fluxuncap5": 513.7743881699428,
    "magap5": 21.974455096230088,
    "sigmagap5": 5.43414211374817,
    "fluxap10": 1175.65266261521,
    "fluxuncap10": 1027.506059351285,
    "magap10": 19.332961527420412,
    "sigmagap10": 1.0954197393390197,
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
