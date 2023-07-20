"""
Python script containing all IO functions.

All opening/writing of fits files should run via this script.
"""

import copy
import logging
import warnings
from pathlib import Path
from typing import Callable

import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning, AstropyWarning

from mirar.data import Image
from mirar.paths import BASE_NAME_KEY, RAW_IMG_KEY, core_fields

logger = logging.getLogger(__name__)


class MissingCoreFieldError(KeyError):
    """Base class for missing core field errors"""


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
    Wrapper function to save an astropy hdu to file

    :param hdu: hdu to save
    :param path: path to save
    :param overwrite: boolean whether to overwrite
    :return: None
    """
    hdu.verify("silentfix+exception")
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


def save_mef_to_path(data_list, header_list, primary_header, path):
    """
    Function to save a MEF image with <data> and <header> to <path>.
    """
    primary_hdu = fits.PrimaryHDU(header=primary_header)
    hdu_list = [primary_hdu]
    assert len(data_list) == len(header_list)
    for ind, data in enumerate(data_list):
        hdu_list.append(fits.ImageHDU(data=data, header=header_list[ind]))

    hdulist = fits.HDUList(hdu_list)

    hdulist.writeto(path, overwrite=True)


def open_fits(path: str | Path) -> tuple[np.ndarray, fits.Header]:
    """
    Function to open a fits file saved to <path>

    :param path: path of fits file
    :return: tuple containing image data and image header
    """
    with fits.open(path, memmap=False) as img:
        hdu = img.pop(0)
        hdu.verify("silentfix+ignore")
        data = hdu.data
        header = hdu.header

    if BASE_NAME_KEY not in header:
        header[BASE_NAME_KEY] = Path(path).name

    if RAW_IMG_KEY not in header.keys():
        header[RAW_IMG_KEY] = path

    return data, header


def open_raw_image(
    path: str | Path,
    open_f: Callable[[str | Path], tuple[np.ndarray, fits.Header]] = open_fits,
) -> Image:
    """
    Function to open a raw image as an Image object

    :param path: path of raw image
    :param open_f: function to open the raw image
    :return: Image object
    """
    data, header = open_f(path)

    new_img = Image(data.astype(np.float64), header)

    check_image_has_core_fields(new_img)

    return new_img


def open_mef_fits(
    path: str | Path,
) -> tuple[fits.Header, list[np.ndarray], list[fits.Header]]:
    """
    Function to open a MEF fits file saved to <path>

    :param path: path of fits file
    :return: tuple containing image data and image header
    """
    split_data, split_headers = [], []
    with fits.open(path, memmap=False) as hdu:
        primary_header = hdu[0].header  # pylint: disable=no-member
        num_ext = len(hdu)
        for ext in range(1, num_ext):
            split_data.append(
                hdu[ext].data.astype(np.float64)
            )  # pylint: disable=no-member
            split_headers.append(hdu[ext].header)  # pylint: disable=no-member

    return primary_header, split_data, split_headers


def combine_mef_extension_file_headers(primary_header, extension_header):
    """
    Function to combine the primary header with an extension header in a MEF frame
    """
    zipped = list(zip(primary_header.values(), primary_header.comments))

    for k in ["XTENSION", "BITPIX"]:
        if k in extension_header.keys():
            del extension_header[k]

    # append primary_header to hdrext
    for count, key in enumerate(list(primary_header.keys())):
        value = zipped[count][0]
        comment = zipped[count][1]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=AstropyWarning)
            if key not in extension_header.keys():
                extension_header.append((key, value, comment))

    return extension_header


def open_mef_image(
    path: str | Path,
    open_f: Callable[
        [str | Path], tuple[fits.Header, list[np.ndarray], list[fits.Header]]
    ] = open_mef_fits,
    extension_key: str | None = None,
) -> list[Image]:
    """
    Function to open a raw image as an Image object

    :param path: path of raw image
    :param open_f: function to open the raw image
    :param extension_key: key to use to number the MEF frames
    :return: Image object
    """

    primary_header, ext_data_list, ext_header_list = open_f(path)

    ext_data_list = [x.astype(np.float64) for x in ext_data_list]
    split_images_list = []

    for ext_num, ext_data in enumerate(ext_data_list):
        ext_header = ext_header_list[ext_num]

        if extension_key is not None:
            extension_num_str = str(ext_header[extension_key])
        else:
            extension_num_str = str(ext_num)

        # append primary_header to hdrext
        new_single_header = combine_mef_extension_file_headers(
            primary_header=primary_header, extension_header=ext_header
        )

        new_single_header[BASE_NAME_KEY] = (
            f"{primary_header[BASE_NAME_KEY].split('.fits')[0]}_"
            f"{extension_num_str}.fits"
        )

        split_images_list.append(
            Image(data=copy.deepcopy(ext_data), header=copy.deepcopy(new_single_header))
        )

    for image in split_images_list:
        check_image_has_core_fields(image)

    return split_images_list


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
                    hdul[-1]._file.tell()  # pylint: disable=protected-access
                    == hdul[-1]._file.size  # pylint: disable=protected-access
                )
        except OSError:
            pass

    return check


def check_image_has_core_fields(img: Image):
    """
    Function to ensure that an image has all the core fields

    :param img: Image object to check
    :return: None
    """
    for key in core_fields:
        if key not in img.keys():
            if BASE_NAME_KEY in img.keys():
                msg = f"({img[BASE_NAME_KEY]}) "

            err = (
                f"New image {msg}is missing the core field {key}. "
                f"Available fields are {list(img.keys())}."
            )
            logger.error(err)
            raise MissingCoreFieldError(err)
