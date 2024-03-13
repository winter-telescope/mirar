"""
Test pipeline visualisation
"""

import logging

from mirar.testing import BaseTestCase
from mirar.utils.pipeline_visualisation import iterate_flowify

logger = logging.getLogger(__name__)


class TestPipelineVisualisation(BaseTestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipeline(self):
        self.logger.info("\n\n Testing pipeline visualisation \n\n")
        iterate_flowify()
