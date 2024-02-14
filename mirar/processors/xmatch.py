"""
Module to cross-match a candidate_table with different catalogs
"""

import logging

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from mirar.catalog.base_catalog import BaseXMatchCatalog
from mirar.data import SourceBatch
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class XMatch(BaseSourceProcessor):
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

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_list in batch:
            candidate_table = source_list.get_data()

            ras = candidate_table["ra"]
            decs = candidate_table["dec"]
            crds = SkyCoord(ras, decs, unit=u.deg)
            query_names = np.array([f"q{x}" for x in np.arange(len(ras))])

            catalog = self.catalog
            logger.debug(f"Querying {catalog.catalog_name} for {len(ras)} sources.")
            query_coords = {
                f"{query_names[ind]}": [ras[ind], decs[ind]] for ind in range(len(ras))
            }
            query_results = catalog.query(query_coords)

            available_projection_keys = []
            for k in catalog.projection.keys():
                if catalog.projection[k] == 1:
                    available_projection_keys += [k]

            # Add placeholder columns for each catalog column
            for key in available_projection_keys:
                for num in range(self.catalog.num_sources):
                    colname = catalog.column_names[key]
                    candidate_table[colname + f"{num + 1}"] = np.array(
                        np.nan,
                        dtype=catalog.column_dtypes[colname],
                    )

            # Add column for number of matches
            nmatch_colname = f"nmtch{self.catalog.abbreviation}"
            candidate_table[nmatch_colname] = 0

            for query_ind, query_name in enumerate(query_names):
                results = query_results[query_name]
                for result_ind, result in enumerate(results):
                    for key in result.keys():
                        colname = catalog.column_names[key] + f"{result_ind + 1}"
                        candidate_table.at[query_ind, colname] = result[key]

                candidate_table.at[query_ind, nmatch_colname] = len(results)

            # Calculate distances between query and result and add to table
            for num in range(self.catalog.num_sources):
                result_ra_colname = self.catalog.ra_column_name + f"{num + 1}"
                result_dec_colname = self.catalog.dec_column_name + f"{num + 1}"
                dist_colname = f"dist{self.catalog.abbreviation}nr{num + 1}"
                candidate_table[dist_colname] = np.array(np.nan, dtype=float)
                crd_nanmask = np.invert(np.isnan(candidate_table[result_ra_colname]))
                result_crds = SkyCoord(
                    ra=candidate_table[result_ra_colname][crd_nanmask],
                    dec=candidate_table[result_dec_colname][crd_nanmask],
                    unit=u.deg,
                )
                candidate_table[dist_colname][crd_nanmask] = (
                    crds[crd_nanmask].separation(result_crds).arcsec
                )

            candidate_table = candidate_table.replace({np.nan: None})

            source_list.set_data(candidate_table)

        return batch
