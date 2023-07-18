"""
Module for loading raw SEDMv2 images and ensuring they have the correct format
"""
import logging
import os
from pathlib import Path

import astropy
import numpy as np

from mirar.data import Image
from mirar.io import open_mef_fits, open_mef_image
from mirar.paths import (
    BASE_NAME_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    __version__,
)

logger = logging.getLogger(__name__)


def clean_science_header(header: astropy.io.fits.Header) -> astropy.io.fits.Header:
    """
    function to modify the primary header of an SEDMv2 science file
    :param header: original primary header of science file
    :return: modified header
    """

    header["OBSTYPE"] = "SCIENCE"
    header["OBSCLASS"] = "science"
    header["TARGET"] = header["OBSTYPE"].lower()

    # positions
    header["RA"] = header["RAD"]
    header["DEC"] = header["DECD"]
    header["TELRA"] = header["TELRAD"]
    header["TELDEC"] = header["TELDECD"]

    header.append(("GAIN", 1.0, "Gain in electrons / ADU"), end=True)
    if header["GAIN"] == 0.0:
        header["GAIN"] = 1.0

    # filters
    header["FILTERID"] = header["FILTER"].split(" ")[1][0]
    header["FILTER"] = header["FILTERID"]

    # keys, IDs ("core fields")
    header[PROC_HISTORY_KEY] = ""
    header["PROCFLAG"] = 0
    header[PROC_FAIL_KEY] = ""
    # header["OBSID"] = 0
    header["PROGID"] = int(3)  # sedmv2's ID
    # header["FIELDID"] = 999999999
    header["COADDS"] = 1
    # header["BZERO"] = 0

    # times
    # orig_dateobs = header["DATE-OBS"]
    # header["DATE-OBS"] = convert_to_UTC(orig_dateobs) # this is the date of entire obs
    # whereas utc in individual headers is for
    # header["UTCTIME"] = header["UTC"]
    # header["TIMEUTC"] = header["UTCTIME"]
    # header["OBSDATE"] = int(header["UTC"].split("_")[0])
    # header["NIGHT"] = int(Time(header["DATE"], format="isot").jd)
    ## - int(Time("2018-01-01", format="iso").jd)  # integer value, night 1, night 2...
    # header["EXPMJD"] = header["OBSDATE"]
    # header["DATE-OBS"] = header["OBSDATE"]

    return header


def clean_cal_header(
    hdr0: astropy.io.fits.Header, hdr1: astropy.io.fits.Header, filepath
) -> tuple[astropy.io.fits.Header, list[astropy.io.fits.Header]]:
    """
    function to modify the primary header of an SEDMv2 calibration file (flat or bias)
    :param hdr0: original primary header of calibration file
    :param hdr1: original secondary header of calibration file
    :return: modified headers
    """

    hdr0["OBSTYPE"] = hdr1["IMGTYPE"].upper()
    hdr0["IMGTYPE"] = hdr1["IMGTYPE"]
    hdr0["OBSCLASS"] = "calibration"
    hdr0["TARGET"] = hdr0[
        "OBSTYPE"
    ].lower()  # should target be OBJECTID? maybe for science imgs?

    hdr0["OBSTYPE"] = hdr0["IMGTYPE"].upper()
    hdr0["OBSCLASS"] = "calibration"
    # flat/bias-specific keys
    if hdr0["IMGTYPE"] == "flat":
        filt = filepath.split("flat_s")[1][0]  # sedm-specific file name structure
        hdr0["FILTERID"] = filt
        hdr0["FILTER"] = f"SDSS {filt}' (Chroma)"
    if hdr0["IMGTYPE"] == "bias":
        hdr0["FILTER"] = "SDSS g"  # arbitrary filter for bias

    hdr0["SOURCE"] = "None"
    hdr0["COADDS"] = 1
    hdr0["CALSTEPS"] = ""
    hdr0["PROCFAIL"] = 1

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
    path: str | Path, skip_first=True
) -> tuple[astropy.io.fits.Header, list[np.array], list[astropy.io.fits.Header]]:
    """
    Load mef image
    """
    header, split_data, split_headers = open_mef_fits(path)

    if "IMGTYPE" in header.keys():
        if skip_first:
            split_data = split_data[1:]
            split_headers = split_headers[1:]
        header = clean_science_header(header)
        # reformat UTC values present in split_headers
        # for hdr in split_headers:
        #    orig = hdr["UTC"]
        #    hdr["UTC"] = convert_to_UTC(orig)
    else:
        header, split_headers = clean_cal_header(header, split_headers[0], path)

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
    Function to load iwinter mef images
    :param path: Path to image
    :return: list of images
    """
    return open_mef_image(path, load_raw_sedmv2_mef)


# def convert_to_UTC(val):
#    """
#    function to reformat the time-related header values in raw sedmv2 files
#    ex:
#        :param val: '20230526_034424.333029'
#        :return: '2023-05-26T03:47:04.310'
#    """
#    val = str(val)
#    yr,mo,day = val[:4],val[4:6],val[6:8]
#    hr,minute,sec = val[9:11],val[11:13],val[13:15]
#    decimal=val.split('.')[1]
#    new_str = yr+'-'+mo+'-'+day+'T'+hr+':'+minute+':'+sec+'.'+decimal
#    t = Time(new_str, format='isot', scale='utc')
#    return t.value
