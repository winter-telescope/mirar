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

logger = logging.getLogger(__name__)

GIT_NONLINEAR_LEVEL = 30000


def load_neowise_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw LT image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)
    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = 1.0

    header["FILTER"] = header["BAND"]

    header["DETCOADD"] = 1
    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = GIT_NONLINEAR_LEVEL * header["DETCOADD"]

    header[OBSCLASS_KEY] = "science"
    header["OBJECT"] = "AT2021blu"
    # Apparently for GIT, the images come tagged correctly.
    header[TARGET_KEY] = header["OBJECT"].lower()
    header["DATE-OBS"] = header["DATE"]
    header["MJD-OBS"] = Time(header["DATE-OBS"]).mjd

    header["JD"] = Time(header["DATE-OBS"]).jd
    if header["BAND"] == 1:
        header["REF_PATH"] = (
            "/Users/viraj/winter_data/neowise/AT2021blu/references"
            "/ICORE_0001_W1_mosaic-int_reference.fits"
        )
    if header["BAND"] == 2:
        header["REF_PATH"] = (
            "/Users/viraj/winter_data/neowise/AT2021blu/references"
            "/ICORE_0001_W2_mosaic-int_reference.fits"
        )
    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    if "FILTERID" not in header.keys():
        header["FILTERID"] = header["BAND"]

    header["EXPTIME"] = 1.0
    header["FID"] = header["FILTERID"]

    if "FIELDID" not in header.keys():
        header["FIELDID"] = 99999
    if "PROGPI" not in header.keys():
        header["PROGPI"] = "Kasliwal"
    if "PROGID" not in header.keys():
        header["PROGID"] = 0

    header["CD1_1"] = header["CDELT1"]
    header["CD1_2"] = 0.0
    header["CD2_1"] = 0.0
    header["CD2_2"] = header["CDELT2"]
    header["ZP"] = header["MAGZP"]
    header[ZP_STD_KEY] = 0.05
    data = data.astype(float)
    data[data == 0.0] = np.nan

    header["OBJRAD"] = 160.643083
    header["OBJDECD"] = +34.437389
    # data[data > 40000] = np.nan
    return data, header


def load_stacked_neowise_image(path: str | Path) -> Image:
    """
    Function to load a raw GIT image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_neowise_fits)
