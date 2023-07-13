"""
Module for querying PS1 catalog
"""
import logging

import numpy as np
from astropy.table import Table

from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog

logger = logging.getLogger(__name__)


class PS1(VizierCatalog):
    """
    PanStarrs 1 catalog
    """

    catalog_vizier_code = "II/349"
    abbreviation = "ps1"

    ra_key = "RAJ2000"
    dec_key = "DEJ2000"

    def filter_catalog(self, table: Table) -> Table:
        logger.debug(f"original ps1 table length: {len(table)}")
        logger.debug("removing ps1 sources with SATURATED flag...")
        sat_flag = 4096  # SATURATED value
        column = table[str(self.filter_name) + "Flags"]
        check = (column & sat_flag) / sat_flag
        # check != 0 means this flag is there
        # check == 0 means this flag is not there
        clean_cat = table[np.where(check == 0)[0]]
        logger.debug(f"found {len(clean_cat)} columns without this flag")
        return clean_cat
