"""
Module for loading raw GIT/LT images and ensuring they have the correct format
"""

import logging
from pathlib import Path

import astropy
import numpy as np
from astropy.time import Time

from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import (
    COADD_KEY,
    GAIN_KEY,
    OBSCLASS_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    SATURATE_KEY,
    TARGET_KEY,
    ZP_STD_KEY,
)

git_filter_dict = {"g": 1, "r": 2, "i": 3, "z": 4, "y": 5}

logger = logging.getLogger(__name__)

GIT_NONLINEAR_LEVEL = 30000


def load_raw_git_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw GIT image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)
    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = 1.0
    header["FILTER"] = header["FILTER"].strip().lower()

    # header[SNCOSMO_KEY] = sncosmo_filters[header["FILTER"].lower()]

    if "COADDS" in header.keys():
        header["DETCOADD"] = header["COADDS"]
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = GIT_NONLINEAR_LEVEL * header["DETCOADD"]

    if header["OBJECT"] in ["acquisition", "pointing", "focus", "none"]:
        header[OBSCLASS_KEY] = header["OBJECT"]
    else:
        header[OBSCLASS_KEY] = "science"

    # Apparently for GIT, the images come tagged correctly.
    header[TARGET_KEY] = header["OBJECT"].lower()
    # header["DATE-OBS"] = header["UTSHUT"]
    header["MJD-OBS"] = Time(header["DATE-OBS"]).mjd

    header["JD"] = Time(header["DATE-OBS"]).jd

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    if "FILTERID" not in header.keys():
        header["FILTERID"] = git_filter_dict[header["FILTER"]]

    header["FID"] = header["FILTERID"]

    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0

    header["ZP"] = header["ZP"]
    header[ZP_STD_KEY] = header["ZP_ERR"]
    data = data.astype(float)
    data[data == 0.0] = np.nan

    data[data > 40000] = np.nan
    return data, header


def load_raw_lt_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw LT image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)
    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = 1.0

    header["FILTER"] = header["FILTER1"][-1].lower()

    # header[SNCOSMO_KEY] = sncosmo_filters[header["FILTER"].lower()]

    header["DETCOADD"] = 1
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = GIT_NONLINEAR_LEVEL * header["DETCOADD"]

    if header["OBJECT"] in ["acquisition", "pointing", "focus", "none"]:
        header[OBSCLASS_KEY] = header["OBJECT"]
    else:
        header[OBSCLASS_KEY] = "science"

    # Apparently for GIT, the images come tagged correctly.
    header[TARGET_KEY] = header["OBJECT"].lower()
    # header["DATE-OBS"] = header["UTSHUT"]
    header["MJD-OBS"] = Time(header["DATE-OBS"]).mjd

    header["JD"] = Time(header["DATE-OBS"]).jd

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    if "FILTERID" not in header.keys():
        header["FILTERID"] = git_filter_dict[header["FILTER"]]

    header["FID"] = header["FILTERID"]

    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0

    # header["ZP"] = header["ZP"]
    # header[ZP_STD_KEY] = header["ZP_ERR"]
    data = data.astype(float)
    data[data == 0.0] = np.nan

    # data[data > 40000] = np.nan
    return data, header


def load_raw_git_image(path: str | Path) -> Image:
    """
    Function to load a raw GIT image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_git_fits)


def load_raw_lt_image(path: str | Path) -> Image:
    """
    Function to load a raw LT image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_lt_fits)
