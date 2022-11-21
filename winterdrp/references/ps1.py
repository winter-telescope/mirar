import logging
from glob import glob
from astropy.io import fits
from winterdrp.references.base_reference_generator import BaseReferenceGenerator
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.wcs import WCS
from winterdrp.data import Image

logger = logging.getLogger(__name__)


class PS1Ref(BaseReferenceGenerator):
    abbreviation = "ps1_ref_lookup"

    def __init__(
            self,
            filter_name: str,
    ):
        super(PS1Ref, self).__init__(filter_name)

    def getimages(self, ra, dec, filters="grizy"):

        """Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
        table = Table.read(url, format='ascii')
        return table

    '''
    def geturl(self, ra, dec, size=240, output_size=None, filters="grizy", format="fits", color=False):

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

        if color and format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.getimages(ra, dec, filters=filters)
        # "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=/rings.v3.skycell/2251/055/rings.v3.skyc
        #  ell.2251.055.stk.g.unconv.fits&format=fits&x=145.767310&y=46.165822&size=3819&wcs=1&imagename=cut
        #  out_rings.v3.skycell.2251.055.stk.g.unconv.fits"
        url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               f"ra={ra}&dec={dec}&size={size}&format={format}")
        if output_size:
            url = url + "&output_size={}".format(output_size)
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[np.argsort(flist)]
        if color:
            if len(table) > 3:
                # pick 3 filters
                table = table[[0, len(table) // 2, len(table) - 1]]
            for i, param in enumerate(["red", "green", "blue"]):
                url = url + "&{}={}".format(param, table['filename'][i])
        else:
            urlbase = url + "&red="
            url = []
            for filename in table['filename']:
                url.append(urlbase + filename)
        return url
    '''

    def geturl(self, ra, dec, size=240, output_size=None, filters="grizy", format="fits", color=False):

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

        if color and format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.getimages(ra, dec, filters=filters)
        # "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=/rings.v3.skycell/2251/055/rings.v3.skyc
        #  ell.2251.055.stk.g.unconv.fits&format=fits&x=145.767310&y=46.165822&size=3819&wcs=1&imagename=cut
        #  out_rings.v3.skycell.2251.055.stk.g.unconv.fits"

        # url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
        #        f"ra={ra}&dec={dec}&size={size}&format={format}")

        urlbase = f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?&red="
        # f"ra={ra}&dec={dec}&size={size}&format={format}")
        url = []
        for filename in table['filename']:
            full_url = urlbase + filename
            full_url += f'&format={format}&x={ra}&y={dec}&size={int(size)}&wcs=1'
            if output_size:
                full_url = full_url + "&output_size={}".format(output_size)

            url.append(full_url)
        return url

    def get_reference(
            self,
            image: Image
    ) -> fits.PrimaryHDU:
        header = image.get_header()

        nx, ny = header['NAXIS1'], header['NAXIS2']
        dx, dy = header['CD1_1'], header['CD2_2']
        img_size_deg = np.max([dx * nx, dy * ny])
        img_size_arcsec = img_size_deg * 3600

        w = WCS(header)
        ra_cent, dec_cent = w.all_pix2world(nx/2, ny/2, 0)

        ps1_img_size_pix = img_size_arcsec / 0.25
        fitsurl = self.geturl(ra_cent, dec_cent, size=ps1_img_size_pix, filters=self.filter_name, format="fits")
        logger.info(fitsurl)
        with fits.open(fitsurl[0]) as hdul:
            refHDU = hdul.copy()[0]
        refHDU.header['GAIN'] = refHDU.header['CELL.GAIN']
        refHDU.header['ZP'] = refHDU.header['FPA.ZP']
        del refHDU.header['HISTORY']
        return refHDU
