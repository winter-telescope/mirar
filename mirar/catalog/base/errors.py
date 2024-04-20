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
