"""
Module for loading raw WINTER images and ensuring they have the correct format
"""
import logging
import os
import warnings
from pathlib import Path

import astropy
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

from mirar.data import Image
from mirar.io import open_mef_fits
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
)
from mirar.pipelines.winter.constants import imgtype_dict, subdets, winter_filters_map
from mirar.pipelines.winter.models import DEFAULT_FIELD, default_program, itid_dict

logger = logging.getLogger(__name__)


def clean_header(header: astropy.io.fits.Header) -> astropy.io.fits.Header:
    """
    Function to clean the header of an image, adding in missing keys and
    correcting values where necessary

    :param header: Header to clean
    :return: Updated header
    """
    header["UTCTIME"] = Time(header["UTCISO"], format="iso").isot

    header["MJD-OBS"] = Time(header["UTCTIME"]).mjd
    header["DATE-OBS"] = Time(header["UTCTIME"]).isot

    header["OBSCLASS"] = ["science", "calibration"][
        header["OBSTYPE"] in ["DARK", "FLAT"]
    ]

    header["EXPTIME"] = np.rint(header["EXPTIME"])

    header["TARGET"] = header["OBSTYPE"].lower()

    if header["TARGNAME"] == "":
        header["TARGNAME"] = f"field_{header['FIELDID']}"

    if (header["FILTERID"] == "dark") & (header["OBSTYPE"] != "BIAS"):
        header["OBSTYPE"] = "DARK"
        header["TARGET"] = "dark"

    header["RA"] = header["RADEG"]
    header["DEC"] = header["DECDEG"]

    obstime = Time(header["UTCTIME"])
    header["EXPID"] = int((obstime.mjd - 59000.0) * 86400.0)  # seconds since 60000 MJD

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = False

    filter_dict = {"J": 1, "H": 2, "Ks": 3}

    if "FILTERID" not in header.keys():
        header["FILTERID"] = filter_dict[header["FILTER"]]
    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0

    if "CTYPE1" not in header:
        header["CTYPE1"] = "RA---TAN"
    if "CTYPE2" not in header:
        header["CTYPE2"] = "DEC--TAN"

    header["FILTER"] = header["FILTERID"]
    if header["FILTER"] == "Hs":
        header["FILTER"] = "H"

    header["FID"] = int(winter_filters_map[header["FILTERID"]])
    logger.debug(f"Obstime is {obstime}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        header["NIGHTDATE"] = obstime.to_datetime().strftime("%Y-%m-%d")

    logger.debug(f"Nightdate is {header['NIGHTDATE']}")

    header["IMGTYPE"] = header["OBSTYPE"]

    if not header["IMGTYPE"] in imgtype_dict:
        header["ITID"] = itid_dict[imgtype_dict["OTHER"]]
    else:
        header["ITID"] = itid_dict[imgtype_dict[header["OBSTYPE"]]]

    header["EXPMJD"] = header["MJD-OBS"]

    for key in ["MOONRA", "MOONDEC", "MOONILLF", "MOONPHAS", "SUNAZ"]:
        if header[key] == "":
            header[key] = 99

    if header["FIELDID"] < 0:
        header["FIELDID"] = DEFAULT_FIELD
    if header["PROGID"] < 0:
        header["PROGID"] = default_program.progid
    # TODO: Get puid from PROGID
    header["PUID"] = header["PROGID"]
    return header


def load_proc_winter_image(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Load proc image
    """
    logger.debug(f"Loading {path}")
    with fits.open(path) as img:
        data = img[0].data  # pylint: disable=no-member
        header = img[0].header  # pylint: disable=no-member
        if "weight" in path:
            header["OBSTYPE"] = "weight"

        header["FILTER"] = header["FILTERID"]

    return data, header


def load_stacked_winter_image(
    path: str | Path,
) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Load proc image
    """
    logger.debug(f"Loading {path}")
    with fits.open(path) as img:
        data = img[0].data  # pylint: disable=no-member
        header = img[0].header  # pylint: disable=no-member
        if "weight" in path:
            header["OBSTYPE"] = "weight"

            header["OBSCLASS"] = "weight"
            header["COADDS"] = 1
            header["TARGET"] = "science"
            header["CALSTEPS"] = ""
            header["PROCFAIL"] = 1
            header["RAWPATH"] = ""
            header["BASENAME"] = os.path.basename(path)
            header["TARGNAME"] = "weight"
        if "UTCTIME" not in header.keys():
            header["UTCTIME"] = "2023-06-14T00:00:00"

    return data, header


def load_winter_mef_image(
    path: str | Path,
) -> tuple[astropy.io.fits.Header, list[np.array], list[astropy.io.fits.Header]]:
    """
    Load mef image.
    """
    header, split_data, split_headers = open_mef_fits(path)

    header = clean_header(header)

    header[BASE_NAME_KEY] = os.path.basename(path)
    header[RAW_IMG_KEY] = path

    # Sometimes there are exptime keys
    for board_header in split_headers:
        if "EXPTIME" in board_header.keys():
            del board_header["EXPTIME"]

    return header, split_data, split_headers


def annotate_winter_subdet_headers(image: Image) -> Image:
    """
    Annotate winter header with information on the subdetector

    :param image: Image to annotate
    :return: Annotated header
    """
    data = image.get_data()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        _, med, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
        image["MEDCOUNT"] = med
        image["STDDEV"] = std
        image["PROCSTATUS"] = 0

    subnx, subny, subnxtot, subnytot = (
        image["SUBNX"],
        image["SUBNY"],
        image["SUBNXTOT"],
        image["SUBNYTOT"],
    )

    mask = (
        (subdets["nx"] == subnx)
        & (subdets["ny"] == subny)
        & (subdets["nxtot"] == subnxtot)
        & (subdets["nytot"] == subnytot)
        & (subdets["boardid"] == image["BOARD_ID"])
    )
    assert np.sum(mask) == 1, (
        f"Subdet not found for nx={subnx}, ny={subny}, "
        f"nxtot={subnxtot}, nytot={subnytot} and boardid={image['BOARD_ID']}"
    )
    image["SUBDETID"] = int(subdets[mask]["subdetid"].iloc[0])
    image["RAWID"] = int(f"{image['EXPID']}_{str(image['SUBDETID']).rjust(2, '0')}")
    image["USTACKID"] = None

    if "DATASEC" in image.keys():
        del image["DATASEC"]

    return image


def get_raw_winter_mask(image: Image) -> np.ndarray:
    """
    Get mask for raw winter image.
    """
    data = image.get_data()
    header = image.header

    mask = np.zeros(data.shape)
    if header["BOARD_ID"] == 0:
        # data[:500, 700:1500] = np.nan
        mask[1075:, :] = 1.0
        mask[:, 1950:] = 1.0
        mask[:20, :] = 1.0

    if header["BOARD_ID"] == 1:
        pass

    if header["BOARD_ID"] == 2:
        mask[1060:, :] = 1.0
        mask[:, 1970:] = 1.0
        mask[:55, :] = 1.0
        mask[:, :20] = 1.0

    if header["BOARD_ID"] == 3:
        mask[1085:, :] = 1.0
        mask[:, 1970:] = 1.0
        mask[:55, :] = 1.0
        mask[:, :20] = 1.0

    if header["BOARD_ID"] == 4:
        # data[610:, :280] = np.nan
        mask[:, 1948:] = 1.0
        mask[:, :61] = 1.0
        mask[:20, :] = 1.0
        mask[1060:, :] = 1.0
        mask[:, 999:1002] = 1.0

    if header["BOARD_ID"] == 5:
        # data[740:, 1270: 1850] = np.nan
        mask[1072:, :] = 1.0
        mask[:, 1940:] = 1.0
        mask[:15, :] = 1.0
        mask[680:, 1180:] = 1.0

    return mask.astype(bool)
