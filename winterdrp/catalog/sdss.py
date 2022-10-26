import logging
import astropy.table
from winterdrp.catalog.base_catalog import VizierCatalog
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from winterdrp.errors import ProcessorError
import randomsdss

logger = logging.getLogger(__name__)


class NotInSDSSError(ProcessorError):
    pass


def in_sdss(
        ra_deg: float,
        dec_deg: float
) -> bool:
    dr16 = randomsdss.DR16()
    return dr16.contains(ra_deg, dec_deg)


class SDSS(VizierCatalog):

    catalog_vizier_code = "V/154"
    abbreviation = "sdss"

    ra_key = "RA_ICRS"
    dec_key = "DE_ICRS"

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        if not in_sdss(ra_deg, dec_deg):
            err = f"Querying for SDSS sources, but the field ({ra_deg}, {dec_deg}) was not observed in SDSS."
            logger.error(err)
            raise NotInSDSSError(err)
