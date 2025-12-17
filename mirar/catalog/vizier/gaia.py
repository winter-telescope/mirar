"""
Module for querying Gaia catalog
"""

import logging

import astropy.table

from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog

logger = logging.getLogger(__name__)

# Disable astroquery warnings
logging.getLogger("astroquery").setLevel(logging.WARNING)


class GaiaVizier(VizierCatalog):
    """
    Gaia DR3 catalog
    """

    catalog_vizier_code = ["I/355/gaiadr3"]

    ra_key = "RA_ICRS"
    dec_key = "DE_ICRS"

    @property
    def abbreviation(self):
        return "gaiadr3"

    @property
    def extra_columns(self) -> list[str]:
        """Code of catalog in Vizier"""
        return ["e_Gmag", "PSS", "Varflag", "e_RPmag", "e_BPmag"]

    def get_mag_key(self) -> str:
        """
        Returns the key for mag in table

        :return: Mag key
        """
        return f"{self.filter_name.upper()}mag"

    def get_column_filters(self) -> dict:
        """
        Get the column filters for the catalog

        :return: Column filters
        """
        return {"PSS": "> 0.99"}

    def get_catalog(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        """
        Query the Gaia catalog for sources around a given position

        :param ra_deg: Right ascension in degrees
        :param dec_deg: Declination in degrees
        :return: Table of sources
        """
        cat = super().get_catalog(ra_deg, dec_deg)

        cat["ra_errdeg"] = cat["e_RA_ICRS"] / 3.6e6
        cat["dec_errdeg"] = cat["e_DE_ICRS"] / 3.6e6

        cat.rename_column("RPmag", "phot_rp_mean_mag")
        cat.rename_column("e_RPmag", "phot_rp_mean_mag_error")
        cat.rename_column("BPmag", "phot_bp_mean_mag")
        cat.rename_column("e_BPmag", "phot_bp_mean_mag_error")

        return cat
