import logging
from glob import glob

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.sdss import SDSS

from winterdrp.data import Image
from winterdrp.references import ReferenceImageError
from winterdrp.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)


class SDSSRef(BaseReferenceGenerator):
    abbreviation = "sdss_ref_lookup"

    def __init__(
        self,
        filter_name: str,
    ):
        super(SDSSRef, self).__init__(filter_name)

    def get_reference(self, image: Image) -> fits.PrimaryHDU:

        header = image.get_header()

        nx, ny = header["NAXIS1"], header["NAXIS2"]

        w = WCS(header)

        ra_cent, dec_cent = w.all_pix2world(nx, ny, 0)

        logger.info(f"Querying SDSS image around {ra_cent},{dec_cent}")

        crd = SkyCoord(ra=ra_cent, dec=dec_cent, unit=(u.deg, u.deg))
        rad = 10
        imgs = []

        while rad < 100:
            imgs = SDSS.get_images(
                crd, radius=rad * u.arcsec, band=self.filter_name.lower()
            )
            if imgs is not None:
                break
            logger.info(
                f"No source found within {rad} arcsec, will try with a larger radius"
            )
            rad += 10

        if len(imgs) == 0:
            err = f"Reference image not found from SDSS"
            logger.error(err)
            raise ReferenceImageError(err)
        else:
            refHDU = imgs[0][0].copy()
            refHDU.header["GAIN"] = 1
            refHDU.header["ZP"] = 2.5 * 9  # Unit of the image is nanomaggie
            del refHDU.header["HISTORY"]
        return refHDU
