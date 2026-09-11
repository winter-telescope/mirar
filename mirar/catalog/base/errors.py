"""
Module for errors in Catalog
"""

from mirar.errors import ProcessorError


class CatalogError(ProcessorError):
    """
    Class for errors in Catalog
    """


class CatalogCacheError(CatalogError):
    """
    Class for errors in CatalogCache
    """


class CatalogQueryError(CatalogError):
    """
    Class for errors querying an external catalog service (e.g. a network
    failure, or a query which unexpectedly returned no usable data)
    """
