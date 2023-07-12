"""
Module for loading raw SEDMv2 images and ensuring they have the correct format
"""
import logging
import os

import astropy
import numpy as np
from astropy.io import fits
from astropy.time import Time

from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_SAVE_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    __version__,
)

logger = logging.getLogger(__name__)


def load_raw_sedmv2_image(path: str) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw SEDMv2 image

    :param path:- path of file
    :return: data and header of image
    """

    with fits.open(path) as data:
        header = data[0].header  # pylint: disable=no-member

        if "PREPTAG" not in header.keys():
            if "IMGTYPE" in header.keys():
                path_new = prepare_science(path)
            else:  # IMGTYPE not in cal; may change w updated cal files
                path_new = prepare_cal(path)
            data = fits.open(path_new)
            header = data[0].header  # pylint: disable=no-member

        header["TARGET"] = header["OBSTYPE"].lower()
        data[0].data = data[0].data * 1.0  # pylint: disable=no-member
        header.append(("GAIN", 1.0, "Gain in electrons / ADU"), end=True)

        # positions
        header["RA"] = header["RAD"]
        header["DEC"] = header["DECD"]
        header["TELRA"] = header["TELRAD"]
        header["TELDEC"] = header["TELDECD"]

        if header["GAIN"] == 0.0:
            header["GAIN"] = 1.0

        # filters
        header["FILTERID"] = header["FILTER"].split(" ")[1][0]
        header["FILTER"] = header["FILTERID"]

        # keys, IDs
        base_name = os.path.basename(path)
        header[BASE_NAME_KEY] = base_name
        header["EXPID"] = int("".join(base_name.split("_")[1:3])[2:])
        header[LATEST_SAVE_KEY] = path
        header[RAW_IMG_KEY] = path
        header[PROC_HISTORY_KEY] = ""
        header["PROCFLAG"] = 0
        header[PROC_FAIL_KEY] = ""
        pipeline_version = __version__
        pipeline_version_padded_str = "".join(
            [x.rjust(2, "0") for x in pipeline_version.split(".")]
        )
        header["PROCID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))
        header["OBSID"] = 0
        header["PROGID"] = int(3)  # sedmv2's ID
        header["FIELDID"] = 999999999
        header["COADDS"] = 1
        header["BZERO"] = 0

        # times
        header["UTCTIME"] = header["UTC"]
        header["TIMEUTC"] = header["UTCTIME"]
        header["OBSDATE"] = int(header["UTC"].split("_")[0])
        header["NIGHT"] = int(Time(header["DATE"], format="isot").jd) - int(
            Time("2018-01-01", format="iso").jd
        )  # integer value, night 1, night 2...
        header["EXPMJD"] = header["OBSDATE"]
        header["DATE-OBS"] = header["OBSDATE"]

    return data[0].data, data[0].header  # pylint: disable=no-member


def prepare_science(filepath: str) -> str:
    """
    Additional steps to get sedmv2 science files into working order

    :param filepath: path of file
    :return: path of revised file
    """

    with fits.open(filepath) as file:
        data = file[0].data  # pylint: disable=no-member
        hdr = file[0].header  # pylint: disable=no-member

        hdr["OBSTYPE"] = "SCIENCE"
        hdr["OBSCLASS"] = "science"
        hdr["PREPTAG"] = 0  # label files that have already run through this function

        # print("changing OBJ coord in header for GW Lib!!!")
        # hdr["OBJRAD"] = 229.9802519789
        # hdr["OBJDECD"] = -25.0067323323

        # discern transient vs. stellar based on filename
        if "ZTF2" in filepath:
            hdr["SOURCE"] = "transient"
        else:
            hdr["SOURCE"] = "stellar"

        # save to file with 1 extension
        fits.writeto(filepath, data, hdr, overwrite=True)  # pylint: disable=no-member
    return filepath


def prepare_cal(filepath: str) -> str:
    """
    Additional steps to get sedmv2 calibration files into working order
    - add required header keywords that are currently missing

    :param filepath: path of file
    :return: path of revised file
    """

    with fits.open(filepath) as cal:
        hdr0, hdr1 = cal[0].header, cal[1].header  # pylint: disable=no-member

        hdr0["OBSTYPE"] = hdr1["IMGTYPE"].upper()
        hdr0["IMGTYPE"] = hdr1["IMGTYPE"]
        hdr0["OBSCLASS"] = "calibration"

        # flat/bias-specific keys
        if hdr0["IMGTYPE"] == "flat":
            filt = filepath.split("flat_s")[1][0]  # sedm-specific file name structure
            hdr0["FILTERID"] = filt
            hdr0["FILTER"] = f"SDSS {filt}' (Chroma)"
        if hdr0["IMGTYPE"] == "bias":
            hdr0["FILTER"] = "SDSS g"

        hdr0["SOURCE"] = "None"

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

        # to label files that have already run through this function
        hdr0["PREPTAG"] = 0
        fits.writeto(
            filepath, cal[1].data, hdr0, overwrite=True  # pylint: disable=no-member
        )
    return filepath
