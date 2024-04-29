"""
Module for querying Gaia catalog
"""

import logging

import astropy.table
from astropy.table import join

from mirar.catalog.base.base_gaia import BaseGaia2Mass
from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog

logger = logging.getLogger(__name__)

# Disable astroquery warnings
logging.getLogger("astroquery").setLevel(logging.WARNING)


class Gaia2MassVizier(BaseGaia2Mass, VizierCatalog):
    """
    Gaia DR3 catalog
    """

    catalog_vizier_code = ["I/355/gaiadr3", "II/246/out"]

    ra_key = "RA_ICRS"
    dec_key = "DE_ICRS"

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
        return {"n2MASS": "= 1", "m2MASS": "= 0"}

    def get_source_table(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        """
        Query the Gaia catalog for sources around a given position

        :param ra_deg: Right ascension in degrees
        :param dec_deg: Declination in degrees
        :return: Table of sources
        """
        tmass_cat = VizierCatalog.get_catalog(self, ra_deg, dec_deg)

        tmass_cat.rename_column("Jmag", "j_m")
        tmass_cat.rename_column("Hmag", "h_m")
        tmass_cat.rename_column("Kmag", "k_m")
        tmass_cat = self.convert_to_ab_mag(tmass_cat)

        # Previously, the magnitude was set in vegamag, but now it is in AB mag`
        tmass_cat["magnitude"] = tmass_cat[f"{self.filter_name.lower()}_m"]

        tmass_cat.rename_column("Qflg", "ph_qual")

        tmass_cat["ra_errdeg"] = tmass_cat["e_RA_ICRS"] / 3.6e6
        tmass_cat["dec_errdeg"] = tmass_cat["e_DE_ICRS"] / 3.6e6

        tmass_cat.rename_column("RPmag", "phot_rp_mean_mag")
        tmass_cat.rename_column("BPmag", "phot_bp_mean_mag")
        tmass_cat.rename_column("e_Jmag", "j_msigcom")
        tmass_cat.rename_column("e_Hmag", "h_msigcom")
        tmass_cat.rename_column("e_Kmag", "k_msigcom")
        return tmass_cat

    def join_query(self, query: dict) -> astropy.table.Table:
        """
        Join the two queries together

        :param query:
        :return:
        """
        return join(query[0], query[1], keys="_2MASS")
