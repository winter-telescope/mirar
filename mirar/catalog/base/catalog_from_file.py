"""
Module for Catalog base class
"""

from pathlib import Path

import astropy.table

from mirar.catalog.base.base_catalog import BaseCatalog
from mirar.utils.ldac_tools import get_table_from_ldac


class CatalogFromFile(BaseCatalog):
    """
    Local catalog from file
    """

    abbreviation = "local"

    def __init__(self, catalog_path: str | Path, *args, **kwargs):
        super().__init__(
            min_mag=0,
            max_mag=99,
            filter_name="None",
            search_radius_arcmin=0,
            *args,
            **kwargs,
        )
        self.catalog_path = catalog_path
        if isinstance(self.catalog_path, str):
            self.catalog_path = Path(self.catalog_path)

    def get_catalog(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        catalog = get_table_from_ldac(self.catalog_path)
        return catalog
