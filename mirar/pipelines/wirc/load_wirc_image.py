"""
Module for loading raw WIRC images and ensuring they have the correct format
"""

import logging
from pathlib import Path

import astropy
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
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
    ZP_KEY,
    ZP_STD_KEY,
)
from mirar.processors.skyportal import SNCOSMO_KEY

wirc_filter_dict = {"J": 1, "H": 2, "K": 3}

logger = logging.getLogger(__name__)

WIRC_NONLINEAR_LEVEL = 30000

sncosmo_filters = {
    "j": "cspjs",
    "h": "csphs",
    "k": "cspk",
}


def load_raw_wirc_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw WIRC image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)

    corrupted = False

    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = 1.2
    header["FILTER"] = header["AFT"].split("__")[0][0]

    header[SNCOSMO_KEY] = sncosmo_filters[header["FILTER"].lower()]

    if "COADDS" in header.keys():
        header["DETCOADD"] = header["COADDS"]
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = WIRC_NONLINEAR_LEVEL * header["DETCOADD"]

    if header["OBJECT"] in ["acquisition", "pointing", "focus", "none", "dark", "flat"]:
        header[OBSCLASS_KEY] = header["OBJECT"]
    else:
        header[OBSCLASS_KEY] = "science"

    # Apparently for WIRC, the images come tagged correctly.
    header[TARGET_KEY] = header["OBJECT"].lower()
    if "MJD-OBS" in header.keys():
        header["DATE-OBS"] = Time(header["MJD-OBS"], format="mjd").isot
    else:
        header["DATE-OBS"] = header["UTSHUT"]
        header["MJD-OBS"] = Time(header["UTSHUT"]).mjd

    header["JD"] = Time(header["DATE-OBS"]).jd

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    try:
        crd = SkyCoord(header["RA"], header["DEC"], unit=(u.hour, u.deg))
        header["CRVAL1"] = crd.ra.deg
        header["CRVAL2"] = crd.dec.deg
    except KeyError:
        header["CRVAL1"] = 0
        header["CRVAL2"] = 0

    if "FILTERID" not in header.keys():
        header["FILTERID"] = wirc_filter_dict[header["FILTER"]]

    header["FID"] = header["FILTERID"]

    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0
    if "ZP" not in header.keys():
        if "TMC_ZP" in header.keys():
            header[ZP_KEY] = float(header["TMC_ZP"])
            header[ZP_STD_KEY] = float(header["TMC_ZPSD"])
        if "ZP_AUTO" in header.keys():
            header[ZP_KEY] = float(header["ZP_AUTO"])
            header[ZP_STD_KEY] = float(header["ZP_AUTO_std"])

    for key in ["TELFOCUS", "RA", "DEC"]:
        if key not in header.keys():
            logger.warning(
                f"No '{key}' entry in header for image {path}. "
                f"Setting as corrupted, will ignore."
            )
            corrupted = True

    header["CORRUPT"] = str(corrupted)

    data = data.astype(float)
    data[data == 0.0] = np.nan
    return data, header


def load_wircpipe_stack_fits(
    path: str | Path,
) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw WIRC image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)

    corrupted = False

    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = 1.2
    header["FILTER"] = header["AFT"].split("__")[0][0]

    header[SNCOSMO_KEY] = sncosmo_filters[header["FILTER"].lower()]

    if "COADDS" in header.keys():
        header["DETCOADD"] = header["COADDS"]
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = WIRC_NONLINEAR_LEVEL * header["DETCOADD"]

    if header["OBJECT"] in ["acquisition", "pointing", "focus", "none", "dark", "flat"]:
        header[OBSCLASS_KEY] = header["OBJECT"]
    else:
        header[OBSCLASS_KEY] = "science"

    # Apparently for WIRC, the images come tagged correctly.
    header[TARGET_KEY] = header["OBJECT"].lower()
    if "MJD-OBS" in header.keys():
        header["DATE-OBS"] = Time(header["MJD-OBS"], format="mjd").isot
    else:
        header["DATE-OBS"] = header["UTSHUT"]
        header["MJD-OBS"] = Time(header["UTSHUT"]).mjd

    header["JD"] = Time(header["DATE-OBS"]).jd

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    if "FILTERID" not in header.keys():
        header["FILTERID"] = wirc_filter_dict[header["FILTER"]]

    header["FID"] = header["FILTERID"]

    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0
    if "ZP" not in header.keys():
        if "TMC_ZP" in header.keys():
            header[ZP_KEY] = float(header["TMC_ZP"])
            header[ZP_STD_KEY] = float(header["TMC_ZPSD"])
        if "ZP_AUTO" in header.keys():
            header[ZP_KEY] = float(header["ZP_AUTO"])
            header[ZP_STD_KEY] = float(header["ZP_AUTO_std"])

    data = data.astype(float)
    # data[data == 0.0] = np.nan
    return data, header


def load_raw_wirc_image(path: str | Path) -> Image:
    """
    Function to load a raw WIRC image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_wirc_fits)


def load_wircpipe_stack_image(path: str | Path) -> Image:
    """
    Function to load a raw WIRC image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_wircpipe_stack_fits)
