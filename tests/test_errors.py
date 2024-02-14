"""
Script for testing the error handling in ..module::mirar.errors
"""

import logging

from mirar.data import Dataset, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)

pipeline = get_pipeline(
    instrument="summer", selected_configurations=["test"], night="20220401"
)

expected_error = {
    "processor_name": "mirar.processors.utils.image_loader",
    "contents": [],
    "known_error_bool": True,
    "non_critical_bool": False,
}


class TestErrors(BaseTestCase):
    """Class for testing errors in ..module::mirar.errors"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing summer pipeline \n\n")

        _, errorstack = pipeline.reduce_images(
            Dataset(ImageBatch()), catch_all_errors=True
        )

        errorstack.summarise_error_stack(verbose=True)

        self.assertEqual(len(errorstack.failed_images), 0)
        self.assertEqual(len(errorstack.noncritical_reports), 0)
        self.assertEqual(len(errorstack.reports), 1)

        err = errorstack.reports[0]

        self.assertTrue(isinstance(err.error, ImageNotFoundError))

        for key, exp in expected_error.items():
            self.assertEqual(exp, getattr(err, key))
