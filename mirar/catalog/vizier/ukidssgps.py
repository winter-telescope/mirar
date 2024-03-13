"""
Module containing UKIDSS GPS Vizier catalog
"""

import logging

from astropy.table import Table

from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog

logger = logging.getLogger(__name__)

TIMEOUT = 30


class UkidssGPS(VizierCatalog):
    """
    UKIDSS GPS catalog
    """

    catalog_vizier_code = "II/316/gps6"
    abbreviation = "ukidssgps"

    ra_key = "RAICRS"
    dec_key = "DEICRS"

    def get_mag_key(self) -> str:
        """
        Returns the key for mag in table

        :return: Mag key
        """
        if self.filter_name == "K":
            return "Kmag1"
        return f"{self.filter_name}mag"

    def get_mag_error_key(self) -> str:
        """
        Returns the key for mag error in table

        :return: Mag error key
        """
        if self.filter_name == "K":
            return "e_Kmag1"
        return f"e_{self.get_mag_key()}"

    def filter_catalog(self, table: Table) -> Table:
        table["ra_errdeg"] = 0.1 / 3600
        table["dec_errdeg"] = 0.1 / 3600
        return table
