"""
Module for PS1 reference image generator
"""
import logging

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS

from mirar.data import Image
from mirar.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)


class PS1Ref(BaseReferenceGenerator):
    """
    PS1 Ref generator
    """

    abbreviation = "ps1_ref_lookup"

    def getimages(self, ra_deg: float, dec_deg: float, filters="grizy") -> Table:
        """
        Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = f"{service}?ra={ra_deg}&dec={dec_deg}&filters={filters}"
        table = Table.read(url, format="ascii")
        return table

    def geturl(
        self,
        ra_deg: float,
        dec_deg: float,
        size: int = 240,
        output_size=None,
        filters="grizy",
        file_format="fits",
        color=False,
    ):
        """Get URL for images in the table

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png" or "fits")
        color = if True, creates a color image (only for jpg or png format).
                Default is return a list of URLs for single-filter grayscale images.
        Returns a string with the URL
        """

        if color and file_format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if file_format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.getimages(ra_deg, dec_deg, filters=filters)

        urlbase = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?&red="

        url = []
        for filename in table["filename"]:
            full_url = urlbase + filename
            full_url += (
                f"&format={file_format}&x={ra_deg}&y={dec_deg}"
                f"&size={int(size)}&wcs=1"
            )
            if output_size:
                full_url += f"&output_size={output_size}"

            url.append(full_url)
        return url

    def get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
        header = image.get_header()

        nx, ny = header["NAXIS1"], header["NAXIS2"]
        dx, dy = header["CD1_1"], header["CD2_2"]
        img_size_deg = np.max([dx * nx, dy * ny])
        img_size_arcsec = img_size_deg * 3600

        wcs = WCS(header)
        ra_cent, dec_cent = wcs.all_pix2world(nx / 2, ny / 2, 0)

        ps1_img_size_pix = img_size_arcsec / 0.25
        fitsurl = self.geturl(
            ra_cent,
            dec_cent,
            size=ps1_img_size_pix,
            filters=self.filter_name,
            file_format="fits",
        )
        logger.debug(fitsurl)

        with fits.open(fitsurl[0]) as hdul:
            ref_hdu = hdul[0].copy()

        ref_hdu.header.rename_keyword("PC001001", "PC1_1")
        ref_hdu.header.rename_keyword("PC001002", "PC1_2")
        ref_hdu.header.rename_keyword("PC002001", "PC2_1")
        ref_hdu.header.rename_keyword("PC002002", "PC2_2")

        ref_hdu.header["GAIN"] = ref_hdu.header["CELL.GAIN"]
        ref_hdu.header["ZP"] = ref_hdu.header["FPA.ZP"]
        del ref_hdu.header["HISTORY"]

        return ref_hdu, None
