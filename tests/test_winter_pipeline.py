"""
Tests for WINTER reduction
"""

import logging
import shutil

from mirar.data import Dataset, ImageBatch
from mirar.paths import get_output_dir
from mirar.pipelines import get_pipeline
from mirar.testing import (
    BaseTestCase,
    load_expected_values,
    save_expected_values,
    update_expected_values,
)

logger = logging.getLogger(__name__)

EXPECTED = load_expected_values("winter")
expected_zp = EXPECTED.get("expected_zp", {})
expected_dataframe_values = EXPECTED.get("expected_dataframe_values", {})
expected_dataframe_ids = EXPECTED.get("expected_dataframe_ids", {})
expected_n_candidates = EXPECTED.get("expected_n_candidates", {})


pipeline = get_pipeline(
    instrument="winter", selected_configurations=["test"], night="20230726"
)

logging.basicConfig(level=logging.DEBUG)


# @unittest.skip(
#     "WFAU is down"
# )
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

    def test_pipeline(self):  # pylint: disable=too-many-branches
        """
        Test winter pipeline
        Returns:

        """
        self.check_expected_values("winter", EXPECTED)
        self.logger.info("\n\n Testing winter pipeline \n\n")

        res, _ = pipeline.reduce_images(Dataset([ImageBatch()]), catch_all_errors=False)

        # Cleanup - delete ouptut dir
        output_dir = get_output_dir(dir_root="winter/20230726")
        shutil.rmtree(output_dir)

        # Expect one dataset, for one different sub-boards
        self.assertEqual(len(res[0]), 1)

        source_table = res[0][0]

        new_candidates_table = source_table.get_data()

        if update_expected_values():
            save_expected_values(
                "winter",
                {
                    "expected_zp": {x: source_table[x] for x in expected_zp},
                    "expected_dataframe_values": {
                        x: list(new_candidates_table[x][:10])
                        for x in expected_dataframe_values
                    },
                    "expected_dataframe_ids": {
                        x: list(new_candidates_table[x][:10])
                        for x in expected_dataframe_ids
                    },
                    "expected_n_candidates": len(new_candidates_table),
                },
            )
            return

        for key, value in expected_zp.items():
            if isinstance(value, float):
                self.assertAlmostEqual(value, source_table[key], places=2)
            elif isinstance(value, int):
                self.assertEqual(value, source_table[key])
            else:
                raise TypeError(
                    f"Type for value ({type(value)} is neither float not int."
                )

        candidates_table = source_table.get_data()

        self.assertEqual(len(candidates_table), expected_n_candidates)
        for key, value in expected_dataframe_values.items():
            for ind, val in enumerate(value):
                if val is None:
                    self.assertIsNone(candidates_table.iloc[ind][key])
                else:
                    self.assertAlmostEqual(
                        candidates_table.iloc[ind][key], val, delta=0.05
                    )

        for key, value in expected_dataframe_ids.items():
            for ind, val in enumerate(value):
                self.assertEqual(candidates_table.iloc[ind][key], val)
