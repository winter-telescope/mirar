"""
catalog contains all the relevant catalogs for use in crossmatching.
The basic structure is contained in :class:`winterdrp.catalog.base_catalog`
"""

from winterdrp.catalog.base_catalog import BaseCatalog
from winterdrp.catalog.gaia import Gaia2Mass
from winterdrp.catalog.ps1 import PS1
from winterdrp.catalog.sdss import SDSS
from winterdrp.catalog.skymapper import SkyMapper
