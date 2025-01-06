"""
Module for loading raw WASP images and ensuring they have the correct format
"""

# pylint: disable=duplicate-code

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
    TIME_KEY,
)
from mirar.pipelines.gmos.config.constants import GMOS_FILTERS

logger = logging.getLogger(__name__)


def load_detrended_gmos_fits(
    path: str | Path,
) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a GMOS image post-DRAGONS processing
    (10.3847/2515-5172/ad0044)

    :param path: path of file
    :return: data and header of image
    """

    primary_header, split_data, split_headers = open_mef_fits(path)

    data = split_data[0]
    header = primary_header
    for key in split_headers[0].keys():
        try:
            header[key] = split_headers[0][key]
        except ValueError:
            pass

    header[RAW_IMG_KEY] = str(path)
    header[BASE_NAME_KEY] = Path(path).name

    if GAIN_KEY not in header.keys():
        header[GAIN_KEY] = header["GAIN"]

    header["FILTER"] = header["FILTER1"].replace(" ", "").split("_")[0].lower()

    assert header["FILTER"] in GMOS_FILTERS, f"Filter {header['FILTER']} not recognised"

    # Set up for forced photometry
    header["OBJRA"] = header["RA"]
    header["OBJDEC"] = header["DEC"]

    del header["RA"]
    del header["DEC"]

    header["DETCOADD"] = header["NCOMBINE"]

    if SATURATE_KEY not in header:
        header[SATURATE_KEY] = header["SATLEVEL"]

    header[OBSCLASS_KEY] = "science"

    header[TARGET_KEY] = header["OBJECT"].lower()

    t = Time(f"{header['DATE-OBS']}T{header['TIME-OBS']}")

    header[TIME_KEY] = t.isot
    header["MJD-OBS"] = t.mjd

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = header["NCOMBINE"]

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = ""

    data = data.astype(float)
    data[data == 0.0] = np.nan

    return data, header


def load_detrended_gmos_image(path: str | Path) -> Image:
    """
    Function to load a detrended GMOS image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_detrended_gmos_fits)
