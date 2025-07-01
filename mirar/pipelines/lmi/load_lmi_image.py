"""
Module for loading raw LMI images and ensuring they have the correct format
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
from mirar.pipelines.lmi.config.constants import (
    LMI_FILTERS,
    LMI_GAIN,
    LMI_NONLINEAR_LEVEL,
)

logger = logging.getLogger(__name__)


def load_raw_lmi_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw LMI image

    :param path: path of file
    :return: data and header of image
    """
    data, header = open_fits(path)
    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = LMI_GAIN

    header["FILTER"] = header["FILTER1"].lower().strip().split("-")[-1]

    assert header["FILTER"] in LMI_FILTERS, f"Filter {header['FILTER']} not recognised"

    # Set up for forced photometry
    ra = Angle(f"{header['OBSRA']} hours").degree
    header["OBJRA"] = ra

    dec = Angle(f"{header['OBSDEC']} degrees").degree
    header["OBJDEC"] = dec

    del header["RA"]
    del header["DEC"]

    header["DETCOADD"] = 1

    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = LMI_NONLINEAR_LEVEL * header["DETCOADD"]

    if header["OBSTYPE"].lower() in ["pointing", "focus", "dark", "bias", "flat"]:
        header[OBSCLASS_KEY] = header["OBSTYPE"].lower()
        header["OBJECT"] = header["OBSTYPE"].lower()
    elif header["OBSTYPE"] == "DOME FLAT":
        header[OBSCLASS_KEY] = "flat"
        header["OBJECT"] = "flat"
    else:
        header[OBSCLASS_KEY] = "science"

    for key in ["focus", "flat", "bias"]:
        if key in Path(path).name.lower():
            header[OBSCLASS_KEY] = key

    header[TARGET_KEY] = header["OBJECT"].lower()

    t = Time(header["DATE-OBS"])

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


def load_raw_lmi_image(path: str | Path) -> Image:
    """
    Function to load a raw LMI image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_lmi_fits)
