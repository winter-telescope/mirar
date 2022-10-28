import logging
from winterdrp.catalog.base_catalog import VizierCatalog
from winterdrp.errors import ProcessorError

logger = logging.getLogger(__name__)


class NotInSkymapperError(ProcessorError):
    pass


def in_skymapper(
        ra_deg: float,
        dec_deg: float
) -> bool:
    return dec_deg < 0.


class SkyMapper(VizierCatalog):

    catalog_vizier_code = "II/358"
    abbreviation = "skymapper"

    ra_key = "RAICRS"
    dec_key = "DEICRS"

    def get_mag_key(self):
        return f"{self.filter_name}PSF"

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        if not in_skymapper(ra_deg, dec_deg):
            err = f"Querying for Skymapper sources, but the field ({ra_deg}, {dec_deg}) was not observed in Skymapper."
            logger.error(err)
            raise NotInSkymapperError(err)
