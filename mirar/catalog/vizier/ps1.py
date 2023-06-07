"""
Module for querying PS1 catalog
"""
import logging

from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog

logger = logging.getLogger(__name__)


class PS1(VizierCatalog):
    """
    PanStarrs 1 catalog
    """

    catalog_vizier_code = "II/349"
    abbreviation = "ps1"

    ra_key = "RAJ2000"
    dec_key = "DEJ2000"
