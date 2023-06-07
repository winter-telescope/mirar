"""
Module for loading raw WIRC images and ensuring they have the correct format
"""
import logging
import os
from pathlib import Path

import astropy
import numpy as np
from astropy.io import fits
from astropy.time import Time

from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
)

logger = logging.getLogger(__name__)


def load_raw_wirc_image(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw WIRC image

    :param path: path of file
    :return: data and header of image
    """
    with fits.open(path) as img:
        # pylint: disable=E1101
        data = img[0].data
        header = img[0].header
        # pylint: enable=E1101
        header["FILTER"] = header["AFT"].split("__")[0]

        if header["OBJECT"] in ["acquisition", "pointing", "focus", "none"]:
            header["OBSTYPE"] = "calibration"

        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]

        header[BASE_NAME_KEY] = os.path.basename(path)
        header[RAW_IMG_KEY] = path
        header["TARGET"] = header["OBJECT"].lower()
        header["UTCTIME"] = header["UTSHUT"]
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
