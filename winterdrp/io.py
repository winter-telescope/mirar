"""
Python script containing all IO functions.

All opening/writing of fits files should run via this script.
"""

import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning


def create_fits(data: np.ndarray, header: fits.Header | None) -> fits.PrimaryHDU:
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


def save_hdu_as_fits(hdu: fits.PrimaryHDU, path: str | Path, overwrite: bool = True):
    """
    Wrapper hunction to save an astropy hdu to file

    :param hdu: hdu to save
    :param path: path to save
    :param overwrite: boolean whether to overwrite
    :return: None
    """
    hdu.writeto(path, overwrite=overwrite)


def save_to_path(
    data: np.ndarray,
    header: fits.Header | None,
    path: str | Path,
    overwrite: bool = True,
):
    """
    Function to save an image with <data> and <header> to <path>.

    :param data: numpy ndarray containing image data
    :param header: astropy Header object
    :param path: output path to save to
    :param overwrite: boolean variable opn whether to overwrite of an
        image exists at <path>. Defaults to True.
    :return: None
    """
    img = create_fits(data, header=header)
    save_hdu_as_fits(hdu=img, path=path, overwrite=overwrite)


def open_fits(path: str | Path) -> tuple[np.ndarray, fits.Header]:
    """
    Function to open a fits file saved to <path>

    :param path: path of fits file
    :return: tuple containing image data and image header
    """
    with fits.open(path) as img:
        hdu = img.pop(0)
        data = hdu.data
        header = hdu.header

    return data, header


def check_file_is_complete(path: str) -> bool:
    """
    Function to check whether a fits file is as large as expected.
    Useful to verify with e.g rsync, where files can be partially transferred

    Disclaimer: I (Robert) do not feel great about having written
    this code block.
    It seems to works though, let's hope no one finds out!
    I will cover my tracks by hiding the astropy warning which
    inspired this block, informing the user that the file
    is not as long as expected

    :param path: path of file to check
    :return: boolean file complete
    """
    check = False
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=AstropyUserWarning)
        try:
            with fits.open(path) as hdul:
                check = (
                    hdul._file.size  # pylint: disable=protected-access
                    == hdul._file.tell()  # pylint: disable=protected-access
                )
        except OSError:
            pass

    return check
