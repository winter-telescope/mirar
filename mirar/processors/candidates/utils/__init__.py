"""
Utils for candidates
"""
import gzip
import io

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from mirar.data import Image
from mirar.processors.candidates.utils.dataframe_writer import DataframeWriter
from mirar.processors.candidates.utils.regions_writer import RegionsWriter


def get_corners_ra_dec_from_header(header: fits.Header) -> list[tuple[float, float]]:
    """
    Function to get corner RA/Dec of the image from the header
    Args:
        header:

    Returns:

    """
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    image_crds = [(0, 0), (nx, 0), (0, ny), (nx, ny)]
    wcs_crds = []
    for image_crd in image_crds:
        ra_deg, dec_deg = WCS(header).all_pix2world([image_crd[0]], [image_crd[1]], 1)
        wcs_crds.append((ra_deg[0], dec_deg[0]))

    return wcs_crds


def get_image_center_wcs_coords(image: Image, origin: int = 0):
    """
    Get RA/Dec of an image
    Args:
        origin: 0 or 1, depending on whether you want it in DS9 coordinates
        image: Image

    Returns:

    """
    header = image.header
    wcs = WCS(header=header)
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    ra_deg, dec_deg = wcs.all_pix2world([nx / 2], [ny / 2], origin)
    return ra_deg[0], dec_deg[0]


def get_xy_from_wcs(
    ra_deg: float | list[float],
    dec_deg: float | list[float],
    header: fits.Header,
    origin: int = 0,
) -> (float, float):
    """
    Get image x-y coordinates of a given ra/dec
    Args:
        ra_deg: RA in decimal degrees
        dec_deg: Dec in decimal degrees
        header: image header, requires WCS keywords
        origin: 0 or 1, depending on whether you want it in DS9 coordinates

    Returns:
        x: x coordinate
        y: y coordinate
    """
    wcs = WCS(header)
    x, y = wcs.all_world2pix(ra_deg, dec_deg, origin)
    return x, y


def makebitims(image: np.array):
    """
    make bit images of the cutouts for the marshal
    Args:
        image: input image cutout

    Returns:
        buf2: a gzipped fits file of the cutout image as
        a BytesIO object
    """
    # open buffer and store image in memory
    buf = io.BytesIO()
    buf2 = io.BytesIO()
    fits.writeto(buf, image)
    with gzip.open(buf2, "wb") as fz:
        fz.write(buf.getvalue())

    return buf2
