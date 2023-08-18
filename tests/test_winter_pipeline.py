"""
Tests for WINTER reduction
"""
import logging
import os
import unittest

from mirar.data import Dataset, ImageBatch
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.176036834716797,
    "ZP_2.0_std": 0.12187111377716064,
    "ZP_2.0_nstars": 721,
    "ZP_3.0": 24.66560935974121,
    "ZP_3.0_std": 0.1159907802939415,
    "ZP_3.0_nstars": 725,
    "ZP_4.0": 24.845422744750977,
    "ZP_4.0_std": 0.1173684149980545,
    "ZP_4.0_nstars": 724,
    "ZP_5.0": 24.897621154785156,
    "ZP_5.0_std": 0.12622451782226562,
    "ZP_5.0_nstars": 729,
    "ZP_6.0": 24.912508010864258,
    "ZP_6.0_std": 0.13024073839187622,
    "ZP_6.0_nstars": 729,
    "ZP_7.0": 24.919605255126953,
    "ZP_7.0_std": 0.13308469951152802,
    "ZP_7.0_nstars": 729,
    "ZP_8.0": 24.924192428588867,
    "ZP_8.0_std": 0.13533064723014832,
    "ZP_8.0_nstars": 730,
    "ZP_AUTO": 24.922786712646484,
    "ZP_AUTO_std": 0.13524207472801208,
    "ZP_AUTO_nstars": 730,
}


pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230710"
)

logging.basicConfig(level=logging.DEBUG)


class TestWinterPipeline(BaseTestCase):
    """
    Module for testing winter pipeline
    """

    def setUp(self):
        """
        Function to set up test
        Returns:

        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        """
        Test winter pipeline
        Returns:

        """
        self.logger.info("\n\n Testing winter pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Expect two datasets, for two different sub-boards
        self.assertEqual(len(res[0]), 1)
        self.assertEqual(len(res[1]), 1)

        headers = [x[0].get_header() for x in res]
        header = [x for x in headers if x["SUBCOORD"] == "0_1"][0]

        # # Uncomment to print new expected ZP dict
        print("New Results WINTER:")
        new_exp = "expected_zp = { \n"
        for header_key in header.keys():
            if header_key in expected_zp.keys():
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
