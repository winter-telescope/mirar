"""
Module for loading raw SEDMv2 images and ensuring they have the correct format
"""
import logging
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from mirar.data import Image
from mirar.io import open_mef_fits, open_mef_image
from mirar.paths import (
    BASE_NAME_KEY,
    GAIN_KEY,
    OBSCLASS_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    TARGET_KEY,
    __version__,
)
from mirar.processors.utils.image_loader import InvalidImage

logger = logging.getLogger(__name__)


def clean_science_header(  # pylint: disable=too-many-branches
    header: fits.Header, split_headers: list[fits.Header], is_mode0: bool
) -> tuple[fits.Header, list[fits.Header]]:
    """
    function to modify the primary header of an SEDMv2 science file
    :param header: original primary header of science file
    :param split_headers: the remaining headers, one for each extension of MEF
    :param is_mode0: True if observed in SEDMv2 observation mode 0
    :return: modified primary header
    """
    if is_mode0:
        # main information is not stored in primary header for mode0
        informative_hdr = split_headers[0]
    else:
        informative_hdr = header

    header[OBSCLASS_KEY] = "science"
    header[TARGET_KEY] = "science"

    # positions
    header["RA"] = informative_hdr["RAD"]
    header["DEC"] = informative_hdr["DECD"]
    header["TELRA"] = informative_hdr["TELRAD"]
    header["TELDEC"] = informative_hdr["TELDECD"]

    # gain
    if GAIN_KEY in informative_hdr:
        if header[GAIN_KEY] == 0.0:
            header[GAIN_KEY] = 1.0
    else:
        header[GAIN_KEY] = 1.0
    for ext in split_headers:
        if GAIN_KEY in ext:
            if ext[GAIN_KEY] == 0.0:
                ext[GAIN_KEY] = 1.0
        else:
            ext[GAIN_KEY] = 1.0

    # filters
    header["FILTERID"] = informative_hdr["FILTER"].split(" ")[1][0]
    header["FILTER"] = header["FILTERID"]
    if is_mode0:
        informative_hdr["FILTER"] = header["FILTER"]
        informative_hdr["FILTERID"] = header["FILTERID"]

    # keys, IDs ("core fields")
    header[PROC_HISTORY_KEY] = ""
    header["PROCFLAG"] = 0
    header[PROC_FAIL_KEY] = ""
    header["PROGID"] = int(3)  # sedmv2's ID
    header["COADDS"] = 1

    if not isinstance(informative_hdr["OBJECTID"], str):
        if isinstance(informative_hdr["QCOMMENT"], str):
            header["OBJECTID"] = informative_hdr["QCOMMENT"].split("_")[0]

    if not "EXPTIME" in informative_hdr:
        header["EXPTIME"] = informative_hdr[
            "EXPOSURE"
        ]  # be weary - exposure is not always reliable

    return header, split_headers


def clean_cal_header(
    hdr0: fits.Header, hdr1: fits.Header, filepath
) -> tuple[fits.Header, list[fits.Header]]:
    """
    function to modify the primary header of an SEDMv2 calibration file (flat or bias)
    :param hdr0: original primary header of calibration file
    :param hdr1: original secondary header of calibration file
    :return: modified headers
    """

    hdr0[OBSCLASS_KEY] = hdr1["IMGTYPE"].lower()
    hdr0["IMGTYPE"] = hdr1["IMGTYPE"]
    hdr0[TARGET_KEY] = hdr0[OBSCLASS_KEY]

    # flat/bias-specific keys
    if hdr0["IMGTYPE"] == "flat":
        filt = filepath.split("flat_s")[1][0]  # sedm-specific file name structure
        hdr0["FILTERID"] = filt
        hdr0["FILTER"] = filt  # f"SDSS {filt}'"  #f"SDSS {filt}' (Chroma)"
    if hdr0["IMGTYPE"] == "bias":
        hdr0["FILTER"] = "g"  # arbitrary filter for bias
        hdr0["FILTERID"] = "g"  # arbitrary filter for bias

    hdr0["SOURCE"] = "None"
    hdr0["COADDS"] = 1
    hdr0["CALSTEPS"] = ""
    hdr0["PROCFAIL"] = 1
    hdr0[GAIN_KEY] = 1.0

    req_headers = [
        "RA",
        "DEC",
        "TELRA",
        "TELDEC",
        "RAD",
        "DECD",
        "TELRAD",
        "TELDECD",
        "UTC",
        "DATE",
    ]
    default_vals = [
        "+00:00:00",
        "+00:00:00",
        "+00:00:00",
        "+00:00:00",
        0.0,
        0.0,
        0.0,
        0.0,
        "20221223_023123.072011",
        "2022-12-23T02:31:23.073",
    ]  # can these be changed? shortened?

    for count, key in enumerate(req_headers):
        hdr0[key] = default_vals[count]

    return hdr0, [hdr1]


def load_raw_sedmv2_mef(
    path: str | Path,
) -> tuple[fits.Header, list[np.array], list[fits.Header]]:
    """
    Load mef image
    """

    sedmv2_ignore_files = [
        "dark",
        "sedm2",
        "speccal",
    ]

    if any(ext in Path(path).name for ext in sedmv2_ignore_files):
        logger.debug(f"Skipping unneeded SEDMv2 file {path}.")
        raise InvalidImage

    header, split_data, split_headers = open_mef_fits(path)

    if "IMGTYPE" in header.keys():  # all modes except mode0
        check_header = header
        skip_first = True
        is_mode0 = False
    else:  # mode0
        check_header = split_headers[0]
        skip_first = False  # there is only 1 extension in mode0
        is_mode0 = True

    if check_header["IMGTYPE"] == "object":
        if skip_first:
            split_data = split_data[1:]
            split_headers = split_headers[1:]
        header, split_headers = clean_science_header(header, split_headers, is_mode0)
    elif check_header["IMGTYPE"] in ["flat", "bias"]:
        header, split_headers = clean_cal_header(header, split_headers[0], path)
    else:
        logger.debug("Unexpected IMGTYPE. Is this a dark?")

    header[BASE_NAME_KEY] = os.path.basename(path)
    header[RAW_IMG_KEY] = path
    header["EXPID"] = int("".join(header[BASE_NAME_KEY].split("_")[1:3])[2:])
    pipeline_version = __version__
    pipeline_version_padded_str = "".join(
        [x.rjust(2, "0") for x in pipeline_version.split(".")]
    )
    header["PROCID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))
    return header, split_data, split_headers


def load_sedmv2_mef_image(
    path: str | Path,
) -> list[Image]:
    """
    Function to load sedmv2 mef images
    :param path: Path to image
    :return: list of images
    """
    return open_mef_image(path, load_raw_sedmv2_mef)
