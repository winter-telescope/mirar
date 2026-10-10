"""
Script for testing the flowchart timings in ..module::mirar.utils.docs
"""

import logging

from mirar.testing import BaseTestCase
from mirar.utils.docs.pipeline_visualisation import format_duration

logger = logging.getLogger(__name__)

expected_durations = {
    0.0: "0.0s",
    9.5: "9.5s",
    59.94: "59.9s",
    59.96: "1m0s",
    60.0: "1m0s",
    61.5: "1m2s",
    599.6: "10m0s",
    3599.4: "59m59s",
    3599.7: "1h0m",
    7265.0: "2h1m",
}


class TestPipelineVisualisation(BaseTestCase):
    """Class for testing flowchart duration formatting"""

    def setUp(self):
        self.logger = logging.getLogger(__name__)

    def test_format_duration(self):
        """
        Test that durations are formatted, and that carries do not leak

        :return: None
        """
        self.logger.info("\n\n Testing duration formatting \n\n")

        for seconds, expected in expected_durations.items():
            self.assertEqual(format_duration(seconds), expected)
