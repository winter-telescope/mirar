"""
Module for loading raw WIRC images and ensuring they have the correct format
"""
import logging
from pathlib import Path

import astropy
import numpy as np
from astropy.time import Time

from mirar.io import open_fits
from mirar.paths import COADD_KEY, PROC_FAIL_KEY, PROC_HISTORY_KEY, SATURATE_KEY

logger = logging.getLogger(__name__)

WIRC_NONLINEAR_LEVEL = 30000


def load_raw_wirc_image(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw WIRC image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)
    header["FILTER"] = header["AFT"].split("__")[0]
    if "COADDS" in header.keys():
        header["DETCOADD"] = header["COADDS"]
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = WIRC_NONLINEAR_LEVEL * header["DETCOADD"]
    if header["OBJECT"] in ["acquisition", "pointing", "focus", "none"]:
        header["OBSTYPE"] = "calibration"

    header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]

    header["TARGET"] = header["OBJECT"].lower()
    if "MJD-OBS" in header.keys():
        header["DATE-OBS"] = Time(header["MJD-OBS"], format="mjd").isot
    else:
        header["DATE-OBS"] = header["UTSHUT"]
        header["MJD-OBS"] = Time(header["UTSHUT"]).mjd
    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    filter_dict = {"J": 1, "H": 2, "Ks": 3}

    if "FILTERID" not in header.keys():
        header["FILTERID"] = filter_dict[header["FILTER"]]
    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0
    if "ZP" not in header.keys():
        if "TMC_ZP" in header.keys():
            header["ZP"] = header["TMC_ZP"]
            header["ZP_std"] = header["TMC_ZPSD"]
        if "ZP_AUTO" in header.keys():
            header["ZP"] = header["ZP_AUTO"]
            header["ZP_std"] = header["ZP_AUTO_std"]
    data = data.astype(float)
    data[data == 0.0] = np.nan
    return data, header
