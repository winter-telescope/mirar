import logging
from glob import glob
from astropy.io import fits
from winterdrp.references.base_reference_generator import BaseReferenceGenerator
import numpy as np
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

logger = logging.getLogger(__name__)


class SDSSRef(BaseReferenceGenerator):
    abbreviation = "sdss_ref_lookup"

    def __init__(
            self,
            filter_name: str,
    ):
        super(SDSSRef, self).__init__(filter_name)

    def get_reference(
            self,
            header: fits.Header
    ) -> fits.PrimaryHDU:
        nx, ny = header['NAXIS1'], header['NAXIS2']

        w = WCS(header)
        ra_cent, dec_cent = w.all_pix2world(nx, ny, 0)

        crd = SkyCoord(ra=ra_cent, dec=dec_cent, unit=(u.deg, u.deg))
        imgs = SDSS.get_images(crd, radius=10 * u.arcsec, band=self.filter_name.lower())
        refHDU = imgs[0][0]
        return refHDU
