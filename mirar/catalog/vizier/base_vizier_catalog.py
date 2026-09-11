"""
Module containing base class for a Vizier catalog
"""

import logging
from abc import ABC

import astropy.table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from mirar.catalog.base.base_catalog import DEFAULT_SNR_THRESHOLD, BaseCatalog
from mirar.catalog.base.errors import CatalogQueryError
from mirar.utils.retry import retry_on_exception

logger = logging.getLogger(__name__)


class VizierError(CatalogQueryError):
    """
    Class for errors in Vizier catalog
    """


class VizierCatalog(BaseCatalog, ABC):
    """
    Base class for a catalog generated using Vizier
    """

    @property
    def catalog_vizier_code(self) -> str | list[str]:
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

    @property
    def extra_columns(self) -> list[str]:
        """Code of catalog in Vizier"""
        return []

    def __init__(
        self,
        *args,
        snr_threshold: float = DEFAULT_SNR_THRESHOLD,
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

    def get_column_filters(self) -> dict:
        """
        Returns the column filters to be applied to the query
        """
        return {}

    @retry_on_exception(exceptions=(VizierError,))
    def get_catalog(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        logger.debug(
            f"Querying {self.abbreviation} catalog around RA {ra_deg:.4f}, "
            f"Dec {dec_deg:.4f} with a radius of {self.search_radius_arcmin:.4f} arcmin"
        )

        viz_cat = Vizier(
            columns=["*"] + self.extra_columns,  # * is only 'all default', not 'all'
            column_filters={
                f"{self.get_mag_key()}": f"{self.min_mag} .. {self.max_mag}",
                f"{self.get_mag_error_key()}": f"<{1.086 / self.snr_threshold:.3f}",
                **self.get_column_filters(),
            },
            row_limit=-1,
            timeout=300,
        )

        try:
            # pylint: disable=no-member
            query = viz_cat.query_region(
                SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg)),
                radius=str(self.search_radius_arcmin) + "m",
                catalog=self.catalog_vizier_code,
                cache=False,
            )
        except Exception as e:
            err = (
                f"Error querying {self.abbreviation} catalog: {e}. "
                "Please check the catalog code and network connection."
            )
            logger.error(err)
            raise VizierError(err) from e

        # An empty result here is ambiguous: it might mean the position is
        # genuinely outside the catalog's footprint (checked below), or it
        # might be a transient/incomplete response from the Vizier service.
        # Either way, retry_on_exception will retry a few times before
        # giving up, since a real "not covered" case will fail identically
        # on every attempt.
        if len(query) == 0:
            self.check_coverage(ra_deg, dec_deg)
            err = f"No matches found in the given radius in {self.abbreviation}"
            logger.error(err)
            raise VizierError(err)

        table = self.join_query(query)

        if len(table) == 0:
            self.check_coverage(ra_deg, dec_deg)
            err = (
                f"Query for {self.abbreviation} returned a table with zero rows "
                f"in the given radius"
            )
            logger.error(err)
            raise VizierError(err)

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
        table["magnitude_err"] = table[self.get_mag_error_key()]
        logger.debug(
            f"{len(table)} matches found in the given radius in {self.abbreviation}"
        )
        table.meta["description"] = ""
        table = self.filter_catalog(table)
        return table

    def join_query(self, query: dict) -> astropy.table.Table:
        """
        Join the query results into a single table

        :param query: Query results
        :return: Table
        """
        return query[0]

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        """
        Perform any available coverage check, to see if catalog covers ra/dec position

        :param ra_deg: Ra
        :param dec_deg: Dec
        :return: None
        """
