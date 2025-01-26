"""
Module for loading raw CFHT images and ensuring they have the correct format
"""

import logging
from pathlib import Path

import astropy
import numpy as np
from astropy.time import Time

from mirar.data import Image
from mirar.io import open_mef_fits, open_raw_image
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    GAIN_KEY,
    OBSCLASS_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    SATURATE_KEY,
    TARGET_KEY,
)

cfht_filter_dict = {"g": 1, "r": 2, "i": 3, "z": 4, "y": 5}

logger = logging.getLogger(__name__)

CFHT_NONLINEARITY_LEVEL = 30000


def load_raw_single_ext_cfht_fits(
    path: str | Path,
) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw CFHT image

    :param path: path of file
    :return: data and header of image
    """
    primary, split_data, split_headers = open_mef_fits(path)
    data = split_data[0]
    header = split_headers[0]
    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = 1.0
    header["FILTER"] = header["FILTER"].strip().lower()

    if "COADDS" in header.keys():
        header["DETCOADD"] = header["COADDS"]
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = CFHT_NONLINEARITY_LEVEL * header["DETCOADD"]

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
        header["FILTERID"] = cfht_filter_dict[header["FILTER"][0]]

    header["FID"] = header["FILTERID"]

    data = data.astype(float)
    # data[data == 0.0] = np.nan
    #
    # data[data > 60000] = np.nan
    header[RAW_IMG_KEY] = str(path)
    header[BASE_NAME_KEY] = path.name
    return data, header


def load_raw_cfht_image(path: str | Path) -> Image:
    """
    Function to load a raw GIT image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_single_ext_cfht_fits)
