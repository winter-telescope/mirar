"""
Base Catalog for crossmatching
"""

from abc import ABC

from mirar.catalog.base.base_catalog import ABCatalog


class BaseXMatchCatalog(ABCatalog, ABC):
    """
    Base Catalog for crossmatching
    """

    @property
    def catalog_name(self):
        """
        Name of catalog
        """
        raise NotImplementedError

    @property
    def projection(self):
        """
        projection for kowalski xmatch
        """
        raise NotImplementedError

    @property
    def column_names(self):
        """
        Name of columns
        """
        raise NotImplementedError

    @property
    def column_dtypes(self):
        """
        dtype of columns
        """
        raise NotImplementedError

    @property
    def ra_column_name(self):
        """
        Name of RA column
        """
        raise NotImplementedError

    @property
    def dec_column_name(self):
        """
        Name of Dec column
        """
        raise NotImplementedError

    def __init__(self, *args, num_sources: int = 1, **kwargs):
        super().__init__(*args, **kwargs)
        self.search_radius_arcsec = self.search_radius_arcmin * 60.0
        self.num_sources = num_sources

    def query(self, coords: dict) -> dict:
        """
        Query coords for result

        :param coords: ra/dec
        :return: crossmatch
        """
        raise NotImplementedError
