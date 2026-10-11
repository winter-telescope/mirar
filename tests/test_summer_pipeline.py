"""
Tests for summer reduction
"""

import logging
import shutil

import pandas as pd

from mirar.data import Dataset, ImageBatch
from mirar.paths import get_output_dir
from mirar.pipelines import get_pipeline
from mirar.testing import BaseTestCase
from mirar.utils.run_stats import RUNTIME_COLUMNS

logger = logging.getLogger(__name__)

expected_zp = {
    "ZP_2.0": 24.376594086837773,
    "ZP_2.0_std": 0.07882569714458094,
    "ZP_2.0_nstars": 30,
    "ZP_3.0": 25.077156103007,
    "ZP_3.0_std": 0.0654024719943491,
    "ZP_3.0_nstars": 30,
    "ZP_4.0": 25.449978212865194,
    "ZP_4.0_std": 0.061531141268491046,
    "ZP_4.0_nstars": 30,
    "ZP_5.0": 25.662833965174357,
    "ZP_5.0_std": 0.0636769631520071,
    "ZP_5.0_nstars": 30,
    "ZP_6.0": 25.781941561253866,
    "ZP_6.0_std": 0.06448472945729802,
    "ZP_6.0_nstars": 30,
    "ZP_7.0": 25.854275071970623,
    "ZP_7.0_std": 0.06542328480566681,
    "ZP_7.0_nstars": 30,
    "ZP_8.0": 25.90200685297648,
    "ZP_8.0_std": 0.06612588202264635,
    "ZP_8.0_nstars": 30,
    "ZP_AUTO": 25.966977791341147,
    "ZP_AUTO_std": 0.07517795091849588,
    "ZP_AUTO_nstars": 30,
}

pipeline = get_pipeline(
    instrument="summer", selected_configurations=["test"], night="20220402"
)


def check_run_stats(test_case: BaseTestCase, stats_df: pd.DataFrame):
    """
    Check the statistics written by a pipeline run

    The runtimes themselves depend on the machine, so only the shape of the table
    and the relations between its columns are checked.

    :param test_case: Test case to assert with
    :param stats_df: Statistics written by the run
    :return: None
    """
    processors = pipeline.all_pipeline_configurations["test"]

    test_case.assertEqual(
        list(stats_df["processor"]), [x.__class__.__name__ for x in processors]
    )
    test_case.assertEqual(list(stats_df["step"]), list(range(1, len(processors) + 1)))
    test_case.assertEqual(set(stats_df["instrument"]), {"summer"})
    test_case.assertEqual(set(stats_df["night"]), {20220402})
    test_case.assertEqual(set(stats_df["configuration"]), {"test"})

    # Nothing should have failed, since the run did not catch errors
    test_case.assertEqual(stats_df["n_errors"].sum(), 0)

    test_case.assertTrue((stats_df[RUNTIME_COLUMNS] >= 0.0).all().all())
    test_case.assertGreater(stats_df["runtime"].sum(), 0.0)

    # A step is normalised by the blocks it handled, which for a loader are its output
    n_blocks = stats_df[["n_input_blocks", "n_output_blocks"]].max(axis=1)
    expected_per_block = (stats_df["runtime"] / n_blocks).where(n_blocks > 0, 0.0)

    for value, expected in zip(stats_df["runtime_per_block"], expected_per_block):
        test_case.assertAlmostEqual(value, expected, delta=0.01)


class TestSummerPipeline(BaseTestCase):
    """
    Module for testing summer pipeline
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
        Test summer pipeline
        Returns:

        """
        self.logger.info("\n\n Testing summer pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Read the run statistics before the output directory holding them is deleted
        stats_df = pd.read_csv(pipeline.get_stats_output_path(), sep="\t")

        # Cleanup - delete non-empty ouptut dir
        output_dir = get_output_dir(dir_root="summer/20220402")
        shutil.rmtree(output_dir)

        check_run_stats(self, stats_df)

        self.assertEqual(len(res[0]), 1)

        header = res[0][0].get_header()

        print("New Results SUMMER:")
        new_exp = "expected_zp = { \n"
        for header_key in header.keys():
            if header_key in expected_zp:
                new_exp += f'    "{header_key}": {header[header_key]}, \n'
        new_exp += "}"
        print(new_exp)

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, header[key], delta=0.01)
            elif isinstance(value, int):
                self.assertEqual(value, header[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )
