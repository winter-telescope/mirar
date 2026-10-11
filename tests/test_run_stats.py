"""
Tests for the run statistics in ..module::mirar.utils.run_stats
"""

import logging
from pathlib import Path

import pandas as pd

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


def make_processor() -> DummyProcessor:
    """
    Build a processor carrying the results of a notional run

    :return: Processor with its latest run recorded
    """
    processor = DummyProcessor()
    processor.latest_runtime = 10.0
    processor.latest_batch_runtimes = [6.0, 9.0]
    processor.latest_n_input_batches = 2
    processor.latest_n_input_blocks = 8
    processor.latest_n_output_batches = 2
    processor.latest_n_output_blocks = 6
    return processor


class TestRunStats(BaseTestCase):
    """
    Class for testing the collection of pipeline run statistics
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)

    def test_latest_stats(self):
        """
        Test that a processor reports its latest run correctly

        :return: None
        """
        self.logger.info("\n\n Testing processor statistics \n\n")

        self.assertEqual(make_processor().latest_stats(), expected_stats)

    def test_n_cpu_capped_by_batches(self):
        """
        Test that the thread count reported is the one actually used

        :return: None
        """
        processor = make_processor()
        processor.latest_n_input_batches = 7

        self.assertEqual(processor.latest_stats()["n_cpu"], processor.max_n_cpu)

    def test_loader_blocks(self):
        """
        Test that a step which creates its own blocks is still normalised

        :return: None
        """
        processor = make_processor()
        processor.latest_n_input_blocks = 0

        self.assertEqual(processor.latest_stats()["runtime_per_block"], 10.0 / 6.0)

    def test_no_blocks(self):
        """
        Test that a step which received nothing does not divide by zero

        :return: None
        """
        processor = DummyProcessor()

        self.assertEqual(processor.latest_stats()["runtime_per_block"], 0.0)
        self.assertEqual(processor.latest_stats()["n_cpu"], 0)

    def test_add_step(self):
        """
        Test that a recorded step carries both the run and the step details

        :return: None
        """
        self.logger.info("\n\n Testing run statistics \n\n")

        run_stats = RunStats(instrument="dummy", night="20230726")
        run_stats.add_step(configuration="test", step=1, processor=make_processor())
        run_stats.add_step(configuration="test", step=2, processor=make_processor())

        stats_df = run_stats.to_dataframe()

        self.assertEqual(len(stats_df), 2)
        self.assertEqual(list(stats_df["step"]), [1, 2])
        self.assertEqual(list(stats_df["configuration"]), ["test", "test"])
        self.assertEqual(set(stats_df["instrument"]), {"dummy"})
        self.assertEqual(set(stats_df["night"]), {"20230726"})
        self.assertEqual(run_stats.total_runtime, 20.0)

        for key, value in expected_stats.items():
            self.assertEqual(stats_df.iloc[0][key], value)

    def test_write_tsv(self):
        """
        Test that the statistics survive a round trip through a tsv

        :return: None
        """
        run_stats = RunStats(instrument="dummy", night="20230726")
        run_stats.add_step(configuration="test", step=1, processor=make_processor())

        output_path = Path(self.temp_dir.name) / "20230726_stats.tsv"
        run_stats.summarise_run_stats(output_path=output_path)

        reloaded = pd.read_csv(output_path, sep="\t")

        self.assertEqual(len(reloaded), 1)
        for key, value in expected_stats.items():
            self.assertEqual(reloaded.iloc[0][key], value)

    def test_empty_run(self):
        """
        Test that a run with no steps still writes a file

        :return: None
        """
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
