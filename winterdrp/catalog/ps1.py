import logging

import astropy.table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.vizier import Vizier

from winterdrp.catalog.base_catalog import VizierCatalog

logger = logging.getLogger(__name__)


class PS1(VizierCatalog):

    catalog_vizier_code = "II/349"
    abbreviation = "ps1"

    ra_key = "RAJ2000"
    dec_key = "DEJ2000"
