"""
Module for catalogs using Vizier
"""

from mirar.catalog.vizier.gaia2mass import Gaia2MassVizier
from mirar.catalog.vizier.ps1 import PS1
from mirar.catalog.vizier.sdss import SDSS, NotInSDSSError, in_sdss
from mirar.catalog.vizier.skymapper import SkyMapper
