import logging
import astropy.table
from astroquery.vizier import Vizier
from winterdrp.catalog.base_catalog import BaseCatalog
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table

logger = logging.getLogger(__name__)


class PS1(BaseCatalog):

    catalog_vizier_code = "II/349"

    def __init__(
            self,
            search_radius_arcmin: float,
            min_mag: float,
            max_mag: float,
            filter_name: str,
            snr_threshold : float=3.0
    ):
        super().__init__(search_radius_arcmin, min_mag, max_mag)
        self.filter_name = filter_name
        self.snr_threshold = snr_threshold

    def get_catalog(
            self,
            ra_deg: float,
            dec_deg: float,
            search_radius_arcmin: float,
            min_mag: float,
            max_mag: float,
    ) -> astropy.table.Table:

        logger.info(
            f'Querying PS1 catalog around RA {ra_deg:.4f}, '
            f'Dec {dec_deg:.4f} with a radius of {search_radius_arcmin:.4f} arcmin'
        )

        v = Vizier(columns=['*'],
                   column_filters={f"{self.filter_name}mag" : f"< {max_mag}",
                                   f"e_{self.filter_name}mag" : "<%.3f" % (1.086 / self.snr_threshold)},
                   row_limit=-1)
        Q = v.query_region(SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg)),
                           radius=str(search_radius_arcmin) + 'm',
                           catalog="II/349", cache=False)

        if len(Q) == 0:
            logger.info('No matches found in the given radius in PS1')
            t = Table()
        else:
            t = Q[0]
            t['RA'] = t['RAJ2000']
            t['DEC'] = t['DEJ2000']
            t['magnitude'] = t[f'{self.filter_name}mag']
            logger.info(f'{len(t)} matches found in the given radius in PS1')

        return t