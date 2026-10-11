"""
Tests for the run statistics in ..module::mirar.utils.run_stats

A real pipeline run is checked in ..module::tests.test_summer_pipeline. These tests
only cover the arithmetic, and the cases a successful run does not reach.
"""

import logging
from pathlib import Path

from mirar.processors.base_processor import BaseProcessor
from mirar.testing import BaseTestCase
from mirar.utils.run_stats import RunStats, get_git_commit

logger = logging.getLogger(__name__)

expected_stats = {
    "processor": "DummyProcessor",
    "n_cpu": 2,
    "runtime": 10.0,
    "batch_runtime_total": 15.0,
    "batch_runtime_max": 9.0,
    "runtime_per_block": 1.25,
    "n_input_batches": 2,
    "n_input_blocks": 8,
    "n_output_batches": 2,
    "n_output_blocks": 6,
    "n_errors": 0,
    "n_noncritical_errors": 0,
}


class DummyProcessor(BaseProcessor):
    """
    Minimal processor, used to check the statistics a processor reports
    """

    base_key = "dummystats"

    max_n_cpu = 4

    def description(self) -> str:
        return "Dummy processor for testing run statistics"

    def _apply(self, batch):
        return batch


class TestRunStats(BaseTestCase):
    """
    Class for testing the collection of pipeline run statistics
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)

    def test_latest_stats(self):
        """
        Test that the statistics are derived correctly from a notional run

        The runtimes of a real run are not known in advance, so this is the only
        place the summed and maximum batch runtimes can be checked against a
        known answer.

        :return: None
        """
        self.logger.info("\n\n Testing processor statistics \n\n")

        processor = DummyProcessor()
        processor.latest_runtime = 10.0
        processor.latest_batch_runtimes = [6.0, 9.0]
        processor.latest_n_input_batches = 2
        processor.latest_n_input_blocks = 8
        processor.latest_n_output_batches = 2
        processor.latest_n_output_blocks = 6

        self.assertEqual(processor.latest_stats(), expected_stats)

    def test_no_blocks(self):
        """
        Test that a step which received nothing does not divide by zero

        :return: None
        """
        processor = DummyProcessor()

        self.assertEqual(processor.latest_stats()["runtime_per_block"], 0.0)
        self.assertEqual(processor.latest_stats()["n_cpu"], 0)

    def test_empty_run(self):
        """
        Test that a run with no steps still writes a file

        :return: None
        """
        self.logger.info("\n\n Testing run statistics \n\n")

        run_stats = RunStats(instrument="dummy", night="20230726")

        output_path = Path(self.temp_dir.name) / "20230726_stats.tsv"
        stats_df = run_stats.summarise_run_stats(output_path=output_path)

        self.assertEqual(len(stats_df), 0)
        self.assertEqual(run_stats.total_runtime, 0)
        self.assertTrue(output_path.exists())

    def test_git_commit(self):
        """
        Test that a git commit can be looked up without raising

        :return: None
        """
        self.assertIsInstance(get_git_commit(), str)
