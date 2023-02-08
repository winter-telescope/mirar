"""
Module for loading raw SEDMv2 images and ensuring they have the correct format
"""
import logging
import os

import astropy
import numpy as np
from astropy.io import fits
from astropy.time import Time

from winterdrp.paths import (
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

    :param path: path of file
    :return: data and header of image
    """
    # print("opening: ", path)
    with fits.open(path) as data:
        header = data[0].header  # pylint: disable=no-member

        if header["IMGTYPE"] == "object":
            header["OBSTYPE"] = "SCIENCE"
        else:
            header["OBSTYPE"] = header["IMGTYPE"].upper()
        header["TARGET"] = header["OBSTYPE"].lower()

        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "SCIENCE"]

        header["RA"] = header["RAD"]
        header["DEC"] = header["DECD"]
        header["TELRA"] = header["TELRAD"]
        header["TELDEC"] = header["TELDECD"]

        # ASSUME THAT CRVAL1 and CRVAL2 are correct based on a-net solution!!!!!!!!
        # header['CRVAL1'] = header['RA']
        # header['CRVAL2'] = header['DEC']

        data[0].data = data[0].data * 1.0  # pylint: disable=no-member

        header["FILTERID"] = header["FILTER"].split(" ")[1][0]
        header["FILTER"] = header["FILTERID"]

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

        header.append(("GAIN", 1.0, "Gain in electrons / ADU"), end=True)

        header["UTCTIME"] = header["UTC"]
        header["TIMEUTC"] = header["UTCTIME"]
        header["OBSDATE"] = int(header["UTC"].split("_")[0])

        # obstime = Time(header["DATE"], format="isot")
        header["NIGHT"] = int(Time(header["DATE"], format="isot").jd) - int(
            Time("2018-01-01", format="iso").jd
        )  # integer value, night 1, night 2...
        header["EXPMJD"] = header["OBSDATE"]

        # IDs
        header["PROGID"] = 0
        header["OBSID"] = 0
        header["SUBPROG"] = "none"
        header["OBSID"] = 0
        header["PROGID"] = "SEDMv2"
        if header["PROGID"] == "SEDMv2":
            header["PROGID"] = 3
        header["PROGID"] = int(header["PROGID"])
        header["FIELDID"] = 999999999
        header["COADDS"] = 1
        header["BZERO"] = 0

    return data[0].data, data[0].header  # pylint: disable=no-member
