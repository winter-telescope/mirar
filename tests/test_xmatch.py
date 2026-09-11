"""
Module for testing the XMatch processor
"""

# pylint: disable=protected-access

import unittest

import pandas as pd

from mirar.catalog.base.base_xmatch_catalog import BaseXMatchCatalog
from mirar.data.source_data import SourceBatch, SourceTable
from mirar.paths import BASE_NAME_KEY, RAW_IMG_KEY
from mirar.processors.xmatch import XMatch


class MockXMatchCatalog(BaseXMatchCatalog):
    """
    Minimal catalog stub for testing XMatch, without needing a real
    external query service.
    """

    catalog_name = "mock_catalog"
    abbreviation = "mock"
    projection = {"_id": 1, "ra": 1, "dec": 1}
    column_names = {"_id": "mockobjectid", "ra": "mockra", "dec": "mockdec"}
    column_dtypes = {"mockobjectid": float, "mockra": float, "mockdec": float}
    ra_column_name = "mockra"
    dec_column_name = "mockdec"

    def __init__(self, match: dict | None):
        super().__init__(search_radius_arcmin=1.0)
        self.match = match

    def query(self, coords: dict) -> dict:
        if self.match is None:
            return {name: [] for name in coords}
        return {name: [self.match] for name in coords}


class TestXMatch(unittest.TestCase):
    """
    Class for testing the XMatch processor
    """

    def make_batch(self) -> SourceBatch:
        """
        Build a minimal single-source batch to cross-match against
        """
        table = pd.DataFrame({"ra": [160.0], "dec": [34.0]})
        metadata = {RAW_IMG_KEY: "raw.fits", BASE_NAME_KEY: "raw.fits"}
        return SourceBatch(SourceTable(table, metadata=metadata))

    def test_large_int_id_returned_as_string(self):
        """
        Some catalog services return large int64-scale ids as JSON strings
        (to avoid precision loss), even though the id column is declared
        as a float. Assigning that string into the float column used to
        be silently accepted (with pandas quietly upcasting the column to
        object dtype); it must still succeed, and the column must remain
        the declared float dtype, not silently degrade.
        """
        catalog = MockXMatchCatalog(
            match={"_id": "924549000121927", "ra": 160.0001, "dec": 34.0001}
        )
        xmatch = XMatch(catalog=catalog)
        batch = xmatch._apply_to_sources(self.make_batch())

        result_table = batch[0].get_data()
        self.assertEqual(result_table.at[0, "mockobjectid1"], 924549000121927.0)
        self.assertEqual(result_table["mockobjectid1"].dtype, float)
        self.assertEqual(result_table.at[0, "nmtchmock"], 1)

    def test_no_match(self):
        """
        A source with no cross-match should get placeholder NaNs and a
        zero match count, not an error.
        """
        catalog = MockXMatchCatalog(match=None)
        xmatch = XMatch(catalog=catalog)
        batch = xmatch._apply_to_sources(self.make_batch())

        result_table = batch[0].get_data()
        self.assertIsNone(result_table.at[0, "mockobjectid1"])
        self.assertEqual(result_table.at[0, "nmtchmock"], 0)


if __name__ == "__main__":
    unittest.main()
