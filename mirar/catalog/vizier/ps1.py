"""
Module for querying PS1 catalog
"""

import logging

import astropy.table
import numpy as np
from astropy.table import Table, join

from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog
from mirar.errors import ProcessorError

logger = logging.getLogger(__name__)


class NotInPS1Error(ProcessorError):
    """Error for source not in PS1 field"""


def in_ps1(dec_deg: float) -> bool:
    """
    Is a given position in PS1?

    :param dec_deg: Declination
    :return: Boolean
    """
    return dec_deg > -30.0


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

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        if not in_ps1(dec_deg):
            err = (
                f"Querying for PS1 sources, but the field "
                f"({ra_deg}, {dec_deg}) was not observed in PS1."
            )
            logger.error(err)
            raise NotInPS1Error(err)


class PS1StarGal(VizierCatalog):
    """
    PanStarrs 1 (PS1) Point Source Catalog (PSC) catalog with Star/Galaxy
    separation by Y. Tachibana & A. A. Miller
    ref: https://iopscience.iop.org/article/10.1088/1538-3873/aae3d9
    """

    catalog_vizier_code = ["II/381/hlsp_ps1_tm", "II/349"]
    abbreviation = "ps1_stargal"

    ra_key = "RAJ2000"
    dec_key = "DEJ2000"

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        if not in_ps1(dec_deg):
            err = (
                f"Querying for PS1_stargal sources, but the field "
                f"({ra_deg}, {dec_deg}) was not observed in PS1_stargal."
            )
            logger.error(err)
            raise NotInPS1Error(err)

    def join_query(self, query: dict) -> astropy.table.Table:
        """
        Join the two queries (PS1 and PS1_TM catalogs) together,
        since PS1_TM only has columns {objid, position, psScore}.

        :param query:
        :return:
        """

        # ps1_tm catalog has worse decimal precision for ra,dec than II/349
        # so we need to remove them before joining
        del query[0]["RAJ2000"]
        del query[0]["DEJ2000"]
        return join(
            query[0],
            query[1],
            keys_left=["objid"],
            keys_right=["objID"],
            join_type="inner",
        )
