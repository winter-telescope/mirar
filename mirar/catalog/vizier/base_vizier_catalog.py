"""
Module containing base class for a Vizier catalog
"""
import logging
from abc import ABC

import astropy.table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.vizier import Vizier

from mirar.catalog.base_catalog import BaseCatalog
from mirar.errors import ProcessorError


class VizierError(ProcessorError):
    """
    Class for errors in Vizier catalog
    """


logger = logging.getLogger(__name__)


class VizierCatalog(BaseCatalog, ABC):
    """
    Base class for a catalog generated using Vizier
    """

    @property
    def catalog_vizier_code(self):
        """Code of catalog in Vizier"""
        raise NotImplementedError()

    @property
    def ra_key(self):
        """Key for RA values"""
        raise NotImplementedError()

    @property
    def dec_key(self):
        """Key for dec values"""
        raise NotImplementedError()

    def __init__(
        self,
        *args,
        snr_threshold: float = 3.0,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.snr_threshold = snr_threshold

    def get_mag_key(self) -> str:
        """
        Returns the key for mag in table

        :return: Mag key
        """
        return f"{self.filter_name}mag"

    def get_mag_error_key(self) -> str:
        """
        Returns the key for mag error in table

        :return: Mag error key
        """
        return f"e_{self.get_mag_key()}"

    def filter_catalog(self, table: astropy.table.Table) -> astropy.table.Table:
        """
        Filters catalog to include a subset of sources, if required
        """
        return table

    def get_catalog(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        logger.debug(
            f"Querying {self.abbreviation} catalog around RA {ra_deg:.4f}, "
            f"Dec {dec_deg:.4f} with a radius of {self.search_radius_arcmin:.4f} arcmin"
        )

        viz_cat = Vizier(
            columns=["*"],
            column_filters={
                f"{self.get_mag_key()}": f"< {self.max_mag}",
                f"{self.get_mag_error_key()}": f"<{1.086 / self.snr_threshold:.3f}",
            },
            row_limit=-1,
        )

        # pylint: disable=no-member
        query = viz_cat.query_region(
            SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg)),
            radius=str(self.search_radius_arcmin) + "m",
            catalog=self.catalog_vizier_code,
            cache=False,
        )

        if len(query) == 0:
            err = f"No matches found in the given radius in {self.abbreviation}"
            logger.error(err)
            self.check_coverage(ra_deg, dec_deg)
            return Table()

        table = query[0]
        logger.debug(f"Table columns are: {table.colnames}")
        if self.get_mag_key() not in table.colnames:
            err = (
                f"Magnitude column {self.get_mag_key()} not found in table."
                f"Available options are : {table.colnames}"
            )
            raise VizierError(err)

        table["ra"] = table[self.ra_key]
        table["dec"] = table[self.dec_key]
        table["magnitude"] = table[self.get_mag_key()]
        logger.debug(
            f"{len(table)} matches found in the given radius in {self.abbreviation}"
        )
        table.meta["description"] = ""
        table = self.filter_catalog(table)
        return table

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        """
        Perform any available coverage check, to see if catalog covers ra/dec position

        :param ra_deg: Ra
        :param dec_deg: Dec
        :return: None
        """
