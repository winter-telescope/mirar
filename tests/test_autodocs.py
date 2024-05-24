"""
Module to check autodocs generation
"""

import logging

from mirar.testing import BaseTestCase
from mirar.utils.docs.auto_config_docs import iterate_rst_generation

logger = logging.getLogger(__name__)


class TestAutodocs(BaseTestCase):
    """
    Class to test automated documentation generation
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipelines(self):
        """
        Test the pipeline configurations

        :return: None
        """
        self.logger.info("\n\n Testing autodocs \n\n")
        iterate_rst_generation()
