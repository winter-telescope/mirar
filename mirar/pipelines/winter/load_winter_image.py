"""
Module for loading raw WINTER images and ensuring they have the correct format
"""
import logging
import os
from pathlib import Path

import astropy
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time

from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
)

logger = logging.getLogger(__name__)


def mask_datasec(data: np.ndarray, header: astropy.io.fits.Header) -> np.ndarray:
    """
    Function to mask the data section of an image
    """
    datasec = header["DATASEC"].replace("[", "").replace("]", "").split(",")
    datasec_xmin = int(datasec[0].split(":")[0])
    datasec_xmax = int(datasec[0].split(":")[1])
    datasec_ymin = int(datasec[1].split(":")[0])
    datasec_ymax = int(datasec[1].split(":")[1])

    data[:, :datasec_xmin] = np.nan
    data[:, datasec_xmax:] = np.nan
    data[:datasec_ymin, :] = np.nan
    data[datasec_ymax:, :] = np.nan
    return data


def load_raw_winter_image(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw WIRC image

    :param path: path of file
    :return: data and header of image
    """
    logger.info(f"Loading {path}")
    with fits.open(path) as img:
        # pylint: disable=E1101
        data = img[0].data
        header = img[0].header

        header["UNIQTYPE"] = f"{header['OBSTYPE']}_{header['BOARD_ID']}"

        basename = os.path.basename(path)
        timestamp = basename.split(".fits")[0].split("_")[1]
        date = timestamp.split("-")[0]
        time = timestamp.split("-")[1]

        if "RADEG" not in header.keys():
            header["RADEG"] = header["RA"]
            header["DECDEG"] = header["DEC"]

        header["UTCTIME"] = (
            f"{date[:4]}-{date[4:6]}-{date[6:]}T" f"{time[:2]}:{time[2:4]}:{time[4:]}"
        )
        logger.info(header["UTCTIME"])
        header["MJD-OBS"] = Time(header["UTCTIME"]).mjd

        header["OBSCLASS"] = ["science", "calibration"][
            header["OBSTYPE"] in ["DARK", "FLAT"]
        ]

        header["EXPTIME"] = np.rint(header["EXPTIME"])
        header[BASE_NAME_KEY] = os.path.basename(path)
        if RAW_IMG_KEY not in header.keys():
            header[RAW_IMG_KEY] = path
        header["TARGET"] = header["OBSTYPE"].lower()

        if header["TARGNAME"] == "":
            header["TARGNAME"] = f"field_{header['FIELDID']}"

        if (header["FILTERID"] == "dark") & (header["OBSTYPE"] != "BIAS"):
            header["OBSTYPE"] = "DARK"
            header["TARGET"] = "dark"

        if ".weight" in path:
            header["OBSTYPE"] = "WEIGHT"
        header["RA"] = header["RADEG"]
        header["DEC"] = header["DECDEG"]

        if COADD_KEY not in header.keys():
            logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
            header[COADD_KEY] = 1

        header[PROC_HISTORY_KEY] = ""
        header[PROC_FAIL_KEY] = ""

        filter_dict = {"J": 1, "H": 2, "Ks": 3}

        if "FILTERID" not in header.keys():
            header["FILTERID"] = filter_dict[header["FILTER"]]
        if "FIELDID" not in header.keys():
            header["FIELDID"] = 99999
        if "PROGPI" not in header.keys():
            header["PROGPI"] = "Kasliwal"
        if "PROGID" not in header.keys():
            header["PROGID"] = 0

        if "CTYPE1" not in header:
            header["CTYPE1"] = "RA---TAN"
        if "CTYPE2" not in header:
            header["CTYPE2"] = "DEC--TAN"

        data = data.astype(float)

        header["FILTER"] = header["FILTERID"]
        if header["FILTER"] == "Hs":
            header["FILTER"] = "H"

        if "DATASEC" in header.keys():
            data = mask_datasec(data, header)
            del header["DATASEC"]
        # TODO: check if this is necessary for short exposures
        #  (Non-linearity/1700 counts)
        data[data > 40000] = np.nan
        if header["BOARD_ID"] == 0:
            # data[:500, 700:1500] = np.nan
            data[1075:, :] = np.nan
            data[:, 1950:] = np.nan
            data[:20, :] = np.nan

        if header["BOARD_ID"] == 1:
            pass

        if header["BOARD_ID"] == 2:
            data[1085:, :] = np.nan
            data[:, 1970:] = np.nan
            data[:55, :] = np.nan
            data[:, :20] = np.nan

        if header["BOARD_ID"] == 3:
            data[1085:, :] = np.nan
            data[:, 1970:] = np.nan
            data[:55, :] = np.nan
            data[:, :20] = np.nan

        if header["BOARD_ID"] == 4:
            # data[610:, :280] = np.nan
            data[:, 1948:] = np.nan
            data[:, :61] = np.nan
            data[:20, :] = np.nan
            data[1060:, :] = np.nan
            data[:, 999:1002] = np.nan

        if header["BOARD_ID"] == 5:
            # data[740:, 1270: 1850] = np.nan
            data[1072:, :] = np.nan
            data[:, 1940:] = np.nan
            data[:15, :] = np.nan
            data[680:, 1180:] = np.nan
            # data[data > 25000] = np.nan

        _, med, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
        header["MEDCOUNT"] = med
        header["STDDEV"] = std
        # if header['RADEG']>300:
        #     header['TARGNAME'] = 'other'

        if "weight" in path:
            header["OBSTYPE"] = "weight"

    return data, header


def load_proc_winter_image(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Load proc image
    """
    logger.info(f"Loading {path}")
    with fits.open(path) as img:
        data = img[0].data
        header = img[0].header
        if "weight" in path:
            header["OBSTYPE"] = "weight"

        header["FILTER"] = header["FILTERID"]

    return data, header


def load_stacked_winter_image(
    path: str | Path,
) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Load proc image
    """
    logger.info(f"Loading {path}")
    with fits.open(path) as img:
        data = img[0].data
        header = img[0].header
        if "weight" in path:
            header["OBSTYPE"] = "weight"

            header["OBSCLASS"] = "weight"
            header["COADDS"] = 1
            header["TARGET"] = "science"
            header["CALSTEPS"] = ""
            header["PROCFAIL"] = 1
            header["RAWPATH"] = ""
            header["BASENAME"] = os.path.basename(path)
            header["TARGNAME"] = "weight"
        if "UTCTIME" not in header.keys():
            header["UTCTIME"] = "2023-06-14T00:00:00"

    return data, header
