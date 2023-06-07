"""
Module to cross-match a candidate_table with different catalogs
"""
import logging

import numpy as np

from mirar.catalog.base_catalog import BaseXMatchCatalog
from mirar.data import SourceBatch
from mirar.processors.base_processor import BaseDataframeProcessor

logger = logging.getLogger(__name__)


class XMatch(BaseDataframeProcessor):
    """
    Class to cross-match a candidate_table to a catalog
    """

    base_key = "XMATCH"

    def __init__(
        self,
        catalog: BaseXMatchCatalog,
    ):
        self.catalog = catalog
        super().__init__()

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_list in batch:
            candidate_table = source_list.get_data()

            ras = candidate_table["ra"]
            decs = candidate_table["dec"]
            query_names = np.array([f"q{x}" for x in np.arange(len(ras))])

            catalog = self.catalog
            logger.info(f"Querying {catalog.catalog_name} for {len(ras)} sources.")
            query_coords = {
                f"{query_names[ind]}": [ras[ind], decs[ind]] for ind in range(len(ras))
            }
            query_results = catalog.query(query_coords)

            available_projection_keys = []
            for k in catalog.projection.keys():
                if catalog.projection[k] == 1:
                    available_projection_keys += [k]

            for key in available_projection_keys:
                for num in range(self.catalog.num_sources):
                    colname = catalog.column_names[key]
                    candidate_table[colname + f"{num + 1}"] = np.array(
                        np.zeros(len(candidate_table)) - 99,
                        dtype=catalog.column_dtypes[colname],
                    )
            nmatch_colname = f"nmtch{self.catalog.abbreviation}"
            candidate_table[nmatch_colname] = 0
            for query_ind, query_name in enumerate(query_names):
                results = query_results[query_name]
                for result_ind, result in enumerate(results):
                    for key in result.keys():
                        colname = catalog.column_names[key] + f"{result_ind + 1}"
                        candidate_table.at[query_ind, colname] = result[key]

                candidate_table.at[query_ind, nmatch_colname] = len(results)

            source_list.set_data(candidate_table)

        return batch
