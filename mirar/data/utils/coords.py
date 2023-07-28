"""
Functions to get image coordinates from WCS
"""
import logging

from astropy.io import fits
from astropy.wcs import WCS

from mirar.data import Image

logger = logging.getLogger(__name__)


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


def write_regions_file(
    regions_path, x_coords, y_coords, system="image", region_radius=5
):
    """
    Function to write a regions file
    Args:
        regions_path: str, path to regions file
        x_coords: list, x-coordinates or RA
        y_coords: list, y-coordinates or Dec
        system: str, image or wcs
        region_radius: float, radius of circle

    Returns:

    """
    logger.debug(f"Writing regions path to {regions_path}")
    with open(f"{regions_path}", "w", encoding="utf8") as regions_f:
        regions_f.write(f"{system}\n")
        for ind, x in enumerate(x_coords):
            regions_f.write(f"CIRCLE({x},{y_coords[ind]},{region_radius})\n")
