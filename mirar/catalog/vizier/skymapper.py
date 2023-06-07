"""
Module containing Skymapper Vizier catalog
"""
import logging

from mirar.catalog.vizier.base_vizier_catalog import VizierCatalog
from mirar.errors import ProcessorError

logger = logging.getLogger(__name__)


class NotInSkymapperError(ProcessorError):
    """Error for source not in Skymapper"""


def in_skymapper(dec_deg: float) -> bool:
    """
    Is a given position in skymapper?

    :param dec_deg: Declination
    :return: Boolean
    """
    return dec_deg < 0.0


class SkyMapper(VizierCatalog):
    """
    Skymapper catalog from Vizier
    """

    catalog_vizier_code = "II/358"
    abbreviation = "skymapper"

    ra_key = "RAICRS"
    dec_key = "DEICRS"

    def get_mag_key(self):
        return f"{self.filter_name}PSF"

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        if not in_skymapper(dec_deg):
            err = (
                f"Querying for Skymapper sources, but the field "
                f"({ra_deg}, {dec_deg}) was not observed in Skymapper."
            )
            logger.error(err)
            raise NotInSkymapperError(err)
