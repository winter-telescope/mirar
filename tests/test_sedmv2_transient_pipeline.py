"""
Module to test SEDMv2 pipeline with "default_transient" configuration
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import get_output_dir
from mirar.pipelines.sedmv2.blocks import process_transient_using_magauto
from mirar.pipelines.sedmv2.load_sedmv2_image import load_sedmv2_mef_image
from mirar.pipelines.sedmv2.sedmv2_pipeline import SEDMv2Pipeline
from mirar.processors.utils import MEFLoader
from mirar.testing import (
    BaseTestCase,
    load_expected_values,
    save_expected_values,
    update_expected_values,
)

logger = logging.getLogger(__name__)

test_data_dir = get_test_data_dir()

EXPECTED = load_expected_values("sedmv2_transient")
expected_zp = EXPECTED.get("expected_zp", {})

test_configuration = [
    MEFLoader(
        input_img_dir=test_data_dir,
        input_sub_dir="raw/",
        load_image=load_sedmv2_mef_image,
    ),
] + process_transient_using_magauto

pipeline = SEDMv2Pipeline(night="20230526", selected_configurations="test_transient")
pipeline.add_configuration(
    configuration_name="test_transient", configuration=test_configuration
)


class TestSEDMv2TransientPipeline(BaseTestCase):
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
        function for testing SEDMv2 transient pipeline
        Returns:

        """
        self.check_expected_values("sedmv2_transient", EXPECTED)
        self.logger.info("\n\n Testing SEDMv2 transient pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Cleanup
        output_dir = get_output_dir("sedmv2/20230526")
        shutil.rmtree(output_dir)

        self.assertEqual(len(res), 1)

        header = res[0][0].get_header()

        if update_expected_values():
            save_expected_values(
                "sedmv2_transient", {"expected_zp": {x: header[x] for x in expected_zp}}
            )
            return

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )
