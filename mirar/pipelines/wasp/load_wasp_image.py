"""
Module for loading raw WASP images and ensuring they have the correct format
"""

# pylint: disable=duplicate-code

import logging
from pathlib import Path

import astropy
import numpy as np
from astropy.coordinates import Angle
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
    TIME_KEY,
)
from mirar.pipelines.wasp.config.constants import (
    WASP_FILTERS,
    WASP_GAIN,
    WASP_NONLINEAR_LEVEL,
)

logger = logging.getLogger(__name__)


def load_raw_wasp_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw LT image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)
    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = WASP_GAIN

    header["FILTER"] = header["FILTER"].lower().replace(" ", "")

    assert header["FILTER"] in WASP_FILTERS, f"Filter {header['FILTER']} not recognised"

    # Set up for forced photometry
    ra = Angle(f"{header['RA']} hours").degree
    header["OBJRA"] = ra

    dec = Angle(f"{header['DEC']} degrees").degree
    header["OBJDEC"] = dec

    del header["RA"]
    del header["DEC"]

    header["DETCOADD"] = 1

    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = WASP_NONLINEAR_LEVEL * header["DETCOADD"]

    if header["OBJECT"] in ["pointing", "focus", "dark", "bias", "flat"]:
        header[OBSCLASS_KEY] = header["OBJECT"]
    else:
        header[OBSCLASS_KEY] = "science"

    for key in ["focus", "flat", "bias"]:
        if key in Path(path).name.lower():
            header[OBSCLASS_KEY] = key

    header[TARGET_KEY] = header["OBJECT"].lower()

    t = Time(header["UTCSTART"])

    header[TIME_KEY] = t.isot
    header["MJD-OBS"] = t.mjd

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    data = data.astype(float)
    data[data == 0.0] = np.nan

    return data, header


def load_raw_wasp_image(path: str | Path) -> Image:
    """
    Function to load a raw WASP image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_wasp_fits)
