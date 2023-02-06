"""
Module for loading raw SEDMv2 images and ensuring they have the correct format
"""
import logging
import os

import astropy
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
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
    with fits.open(path) as data:
        header = data[0].header

        # science / flat / bias / etc...
        if header["IMGTYPE"] == "object":
            header["OBSTYPE"] = "SCIENCE"
        else:
            header["OBSTYPE"] = header["IMGTYPE"].upper()
        header["TARGET"] = header["OBSTYPE"].lower()

        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "SCIENCE"]

        # coordinates
        ## change to deg units
        crd = SkyCoord(
            ra=header["RA"], dec=header["DEC"], unit=(u.hourangle, u.deg)
        )  # HOURANGLE!!!!!
        header["RA"] = crd.ra.deg
        header["DEC"] = crd.dec.deg

        tel_crd = SkyCoord(
            ra=header["TELRA"], dec=header["TELDEC"], unit=(u.hourangle, u.deg)
        )  # HOURANGLE!
        header["TELRA"] = tel_crd.ra.deg
        header["TELDEC"] = tel_crd.dec.deg

        # ASSUME THAT CRVAL1 and CRVAL2 are correct based on a-net solution!!!!!!!!
        # header['CRVAL1'] = header['RA']
        # header['CRVAL2'] = header['DEC']

        data[0].data = data[0].data * 1.0

        # filter
        header["FILTERID"] = header["FILTER"].split(" ")[1][
            0
        ]  # overwrites numerical filterid
        header["FILTER"] = header["FILTERID"]

        # keys...
        header[LATEST_SAVE_KEY] = path
        header[RAW_IMG_KEY] = path
        header[PROC_HISTORY_KEY] = ""
        header["PROCFLAG"] = 0
        header[PROC_FAIL_KEY] = ""

        base_name = os.path.basename(path)
        header[BASE_NAME_KEY] = base_name
        header["EXPID"] = int("".join(base_name.split("_")[1:3])[2:])

        pipeline_version = __version__
        pipeline_version_padded_str = "".join(
            [x.rjust(2, "0") for x in pipeline_version.split(".")]
        )
        header["PROCID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))

        header.append(("GAIN", 1.0, "Gain in electrons / ADU"), end=True)

        # times
        header["UTCTIME"] = header["UTC"]
        header["TIMEUTC"] = header["UTCTIME"]
        header["OBSDATE"] = int(header["UTC"].split("_")[0])

        obstime = Time(header["DATE"], format="isot")
        header["NIGHT"] = int(obstime.jd) - int(
            Time("2018-01-01", format="iso").jd
        )  # integer value, night 1, night 2...
        header["EXPMJD"] = header["OBSDATE"]
        header["SHUTOPEN"] = obstime.jd
        header["SHUTCLSD"] = obstime.jd

        # IDs
        default_id = 0

        for key in ["PROGID", "OBSID"]:
            if key not in header.keys():
                header[key] = default_id
            else:
                try:
                    header[key] = int(header[key])
                except ValueError:
                    header[key] = default_id

        if "SUBPROG" not in header.keys():
            header["SUBPROG"] = "none"

        if "OBSID" not in header.keys():
            logger.warning(f"No {key} found in header of {path}")
            header["OBSID"] = default_id
        else:
            try:
                header["OBSID"] = int(header["OBSID"])
            except ValueError:
                header["OBSID"] = default_id

        header["PROGID"] = "SEDMv2"
        if header["PROGID"] == "SEDMv2":
            header["PROGID"] = 3
        header["PROGID"] = int(header["PROGID"])  # might cause error

        # itid_dict = {
        #    "SCIENCE": 1,
        #    "BIAS": 2,
        #    "FLAT": 2,
        #    "DARK": 2,
        #    "FOCUS": 3,
        #    "POINTING": 4,
        #    "OTHER": 5,
        # }  # may be unnecessary

        # if not header["OBSTYPE"] in itid_dict.keys():
        #    header["ITID"] = 5
        # else:
        #    header["ITID"] = itid_dict[header["OBSTYPE"]]

        header["FIELDID"] = 999999999

        # others
        if "COADDS" not in header.keys():
            header["COADDS"] = 1
        header["BZERO"] = 0

        # sunmoon_keywords = [
        #    "MOONRA",
        #    "MOONDEC",
        #    "MOONILLF",
        #    "MOONPHAS",
        #    "MOONALT",
        #    "SUNAZ",
        #    "SUNALT",
        # ]

        # for key in sunmoon_keywords:
        #    val = 0
        #    if key in header.keys():
        #        if header[key] not in [""]:
        #            val = header[key]
        #    header[key] = val

    return data[0].data, data[0].header
