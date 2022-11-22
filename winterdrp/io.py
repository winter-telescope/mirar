"""
Python script containing all IO functions.

All opening/writing of fits files should run via this script.
"""

from astropy.io import fits
import numpy as np
import astropy.io.fits
from pathlib import Path


def create_fits(
        data: np.ndarray,
        header: astropy.io.fits.Header | None
) -> astropy.io.fits.PrimaryHDU:
    """
    Return an astropy PrimaryHDU object created with <data> and <header>

    :param data: numpy ndarray containing image data
    :param header: astropy Header object
    :return: astropy PrimaryHDU object containing the image data and header
    """
    proc_hdu = fits.PrimaryHDU(data)
    if header is not None:
        proc_hdu.header = header
    return proc_hdu


def save_to_path(
        data: np.ndarray,
        header: astropy.io.fits.Header | None,
        path: str | Path,
        overwrite: bool = True
):
    """
    Function to save an image with <data> and <header> to <path>.

    :param data: numpy ndarray containing image data
    :param header: astropy Header object
    :param path: output path to save to
    :param overwrite: boolean variable opn whether to overwrite of an image exists at <path>. Defaults to True.
    :return: None
    """
    img = create_fits(data, header=header)
    img.writeto(path, overwrite=overwrite)


def open_fits(
        path: str | Path
) -> tuple[np.ndarray, astropy.io.fits.Header]:
    """
    Function to open a fits file saved to <path>

    :param path: path of fits file
    :return: tuple containing image data and image header
    """
    with fits.open(path) as img:
        data = img[0].data
        header = img[0].header

    return data, header
