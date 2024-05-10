"""
catalog contains all the relevant catalogs for use in crossmatching.
The basic structure is contained in :class:`mirar.catalog.base_catalog`
"""

from mirar.catalog.base.base_catalog import BaseCatalog
from mirar.catalog.base.catalog_from_file import CatalogFromFile
from mirar.catalog.multibackend import Gaia2Mass
from mirar.catalog.vizier.ps1 import PS1
