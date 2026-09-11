"""
Module for generating reference images from SDSS
"""

import logging

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.sdss import SDSS

from mirar.data import Image
from mirar.references.base_reference_generator import BaseReferenceGenerator
from mirar.references.errors import ReferenceImageError
from mirar.utils.retry import retry_on_exception

logger = logging.getLogger(__name__)


class SDSSRef(BaseReferenceGenerator):
    """
    SDSS ref generator
    """

    abbreviation = "sdss_ref_lookup"

    @retry_on_exception(exceptions=(ReferenceImageError,))
    def _get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
        header = image.get_header()

        nx, ny = header["NAXIS1"], header["NAXIS2"]

        wcs = WCS(header)

        ra_cent, dec_cent = wcs.all_pix2world(nx / 2, ny / 2, 0)

        logger.debug(f"Querying SDSS image around {ra_cent},{dec_cent}")

        crd = SkyCoord(ra=ra_cent, dec=dec_cent, unit=(u.deg, u.deg))
        rad = 10
        imgs = None

        while rad < 100:
            try:
                imgs = SDSS.get_images(
                    coordinates=crd,
                    radius=rad * u.arcsec,
                    band=self.filter_name.lower(),
                )
            except Exception as e:
                err = f"Error querying SDSS for reference image: {e}."
                logger.error(err)
                raise ReferenceImageError(err) from e
            if imgs is not None:
                break
            logger.debug(
                f"No source found within {rad} arcsec, will try with a larger radius"
            )
            rad += 10

        # imgs may be None (radius exhausted without a hit) or an empty
        # sequence; both mean no reference image was found.
        if not imgs:
            err = f"Reference image not found from SDSS for {self.filter_name.lower()}"
            logger.error(err)
            raise ReferenceImageError(err)

        ref_hdu = imgs[0][0].copy()
        ref_hdu.header["GAIN"] = 1
        ref_hdu.header["ZP"] = 2.5 * 9  # Unit of the image is nanomaggie
        del ref_hdu.header["HISTORY"]

        return ref_hdu, None
