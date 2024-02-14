"""
Module for loading raw WINTER images and ensuring they have the correct format
"""

import logging
import os
import warnings
from pathlib import Path

import astropy
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

from mirar.data import Image, ImageBatch
from mirar.io import (
    open_fits,
    open_mef_fits,
    open_mef_image,
    open_raw_image,
    tag_mef_extension_file_headers,
)
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
    core_fields,
)
from mirar.pipelines.winter.constants import (
    imgtype_dict,
    palomar_observer,
    subdets,
    winter_filters_map,
)
from mirar.pipelines.winter.models import DEFAULT_FIELD, default_program, itid_dict
from mirar.processors.skyportal import SNCOSMO_KEY

logger = logging.getLogger(__name__)

sncosmo_filters = {
    "y": "desy",
    "j": "2massj",
    "h": "2massh",
}


def clean_header(header: fits.Header) -> fits.Header:
    """
    Function to clean the header of an image, adding in missing keys and
    correcting values where necessary

    :param header: Header to clean
    :return: Updated header
    """
    header[GAIN_KEY] = 1.0
    header[SATURATE_KEY] = 40000.0

    header["UTCTIME"] = Time(header["UTCISO"], format="iso").isot

    date_t = Time(header["UTCTIME"])

    header["MJD-OBS"] = date_t.mjd
    header["JD"] = date_t.jd
    header["DATE-OBS"] = date_t.isot

    header[OBSCLASS_KEY] = header["OBSTYPE"].lower()

    # Check to ensure that biases and darks are tagged appropriately
    if header["EXPTIME"] == 0.0:
        header[OBSCLASS_KEY] = "bias"

    if header["FILTERID"] == "dark":
        if header[OBSCLASS_KEY] not in ["dark", "bias", "test", "science"]:
            header[OBSCLASS_KEY] = "test"
        else:
            header[OBSCLASS_KEY] = "dark"

    # Discard pre-sunset, post-sunset darks
    if header[OBSCLASS_KEY] == "dark":
        sun_pos = palomar_observer.sun_altaz(Time(header["UTCTIME"]))
        if sun_pos.alt.to_value("deg") > -20.0:
            header[OBSCLASS_KEY] = "test"

    # Sometimes darks come with wrong fieldids
    if header[OBSCLASS_KEY] == "dark":
        header["FIELDID"] = DEFAULT_FIELD

    header["EXPTIME"] = np.rint(header["EXPTIME"])

    # Set up the target name

    target = f"field_{header['FIELDID']}"
    if ("SCHDNAME" in header.keys()) & ("OBHISTID" in header.keys()):
        if header["SCHDNAME"] != "":
            target = f"{header['SCHDNAME']}_{header['OBHISTID']}"
    elif TARGET_KEY in header.keys():
        if header[TARGET_KEY] != "":
            target = header[TARGET_KEY]
    # If the observation is dark/bias/focus/pointing/flat, enforce TARGET_KEY is also
    # dark/bias as the TARGET_KEY is used for CalHunter. Currently they may not come
    # with the correct TARGET_KEY.
    if header[OBSCLASS_KEY].lower() in [
        "dark",
        "bias",
        "focus",
        "pointing",
        "flat",
        "test",
    ]:
        target = header[OBSCLASS_KEY].lower()

    header[TARGET_KEY] = target

    if "TARGNAME" in header.keys():
        if header["TARGNAME"] == "":
            header["TARGNAME"] = None

    header["RA"] = header["RADEG"]
    header["DEC"] = header["DECDEG"]

    header["EXPID"] = int((date_t.mjd - 59000.0) * 86400.0)  # seconds since 60000 MJD

    if COADD_KEY not in header.keys():
        logger.debug(f"No {COADD_KEY} entry. Setting coadds to 1.")
        header[COADD_KEY] = 1

    header[PROC_HISTORY_KEY] = ""
    header[PROC_FAIL_KEY] = False

    # Make sure filter is a keyword that the pipeline recognizes
    filter_dict = {"J": 1, "H": 2, "Ks": 3}
    if "FILTERID" not in header.keys():
        header["FILTERID"] = filter_dict[header["FILTER"]]
    header["FILTER"] = header["FILTERID"]
    if header["FILTER"] == "Hs":
        header["FILTER"] = "H"
    if header["FILTERID"] in winter_filters_map:
        header["FID"] = int(winter_filters_map[header["FILTERID"]])
    else:
        header["FID"] = -99

    # Set default values if field details seem incorrect
    if "FIELDID" not in header.keys():
        header["FIELDID"] = DEFAULT_FIELD
    if header["FIELDID"] < 0:
        header["FIELDID"] = DEFAULT_FIELD

    # Set default values if program is not correct
    if "PROGPI" not in header.keys():
        header["PROGPI"] = default_program.pi_name
    if "PROGID" not in header.keys():
        header["PROGID"] = default_program.progid
    # If PROGNAME is not present or is empty, set it to default here.
    # Otherwise, it gets set to default in the insert_entry for exposures.
    if "PROGNAME" not in header:
        header["PROGNAME"] = default_program.progname
    if header["PROGNAME"] == "":
        header["PROGNAME"] = default_program.progname

    if len(header["PROGNAME"]) != 8:
        logger.warning(
            f"PROGNAME {header['PROGNAME']} is not 8 characters long. "
            f"Replacing with default."
        )
        header["PROGNAME"] = default_program.progname

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        header["NIGHTDATE"] = date_t.to_datetime().strftime("%Y-%m-%d")

    header["IMGTYPE"] = header[OBSCLASS_KEY]

    if not header["IMGTYPE"] in imgtype_dict:
        header["ITID"] = itid_dict[imgtype_dict["other"]]
    else:
        header["ITID"] = itid_dict[imgtype_dict[header[OBSCLASS_KEY]]]

    header["EXPMJD"] = header["MJD-OBS"]
    if header["FILTER"].lower() in ["y", "j", "h"]:
        header[SNCOSMO_KEY] = sncosmo_filters[header["FILTER"].lower()]

    header["DITHGRP"] = int(header["DITHNUM"] <= 5)
    if "GAINCOLT" not in header.keys():
        header["GAINCOLT"] = "[]"
    if "GAINCOLB" not in header.keys():
        header["GAINCOLB"] = "[]"
    if "GAINROW" not in header.keys():
        header["GAINROW"] = "[]"
    return header


def load_winter_stack(
    path: str | Path,
) -> Image:
    """
    Load proc image

    :param path: Path to image
    :return: data and header
    """
    logger.debug(f"Loading {path}")
    data, header = open_fits(path)

    dirname = path.split("/winter/")[0] + "/winter/"
    wghtpath = header["WGHTPATH"]
    weight_pathname = wghtpath.split("/winter/")[-1]
    new_weightpath = Path(dirname) / weight_pathname
    header["WGHTPATH"] = new_weightpath.as_posix()
    header["SAVEPATH"] = path

    if SNCOSMO_KEY not in header.keys():
        if header["FILTER"].lower() in ["y", "j", "h"]:
            header[SNCOSMO_KEY] = sncosmo_filters[header["FILTER"].lower()]
    return Image(data=data, header=header)


def load_stacked_winter_image(
    path: str | Path,
) -> tuple[np.array, fits.Header]:
    """
    Load proc image

    :param path: Path to image
    :return: data and header
    """
    logger.debug(f"Loading {path}")
    data, header = open_fits(path)

    if "weight" in path:
        header[OBSCLASS_KEY] = "weight"

        header["COADDS"] = 1
        header["CALSTEPS"] = ""
        header["PROCFAIL"] = 1
        header["RAWPATH"] = ""
        header["BASENAME"] = os.path.basename(path)
        header[TARGET_KEY] = "weight"
    if "UTCTIME" not in header.keys():
        header["UTCTIME"] = "2023-06-14T00:00:00"

    return data, header


def load_test_winter_image(
    path: str | Path,
) -> Image:
    """
    Load test WINTER image

    :param path: Path to image
    :return: Image object
    """
    image = open_raw_image(path)
    header = clean_header(image.header)

    image.set_header(header)
    return image


def load_raw_winter_mef(
    path: str | Path,
) -> tuple[astropy.io.fits.Header, list[np.array], list[astropy.io.fits.Header]]:
    """
    Load mef image.

    :param path: Path to image
    :return: Primary header, list of data arrays, list of headers
    """
    primary_header, split_data, split_headers = open_mef_fits(path)

    try:
        primary_header = clean_header(primary_header)
    except KeyError:
        logger.warning(f"Could not clean header for {path}, missing keywords")
        primary_header[OBSCLASS_KEY] = "test"
        primary_header["COADDS"] = 1
        primary_header["CALSTEPS"] = ""
        primary_header["PROCFAIL"] = 1
        primary_header[TARGET_KEY] = "test"
        primary_header["UTCTIME"] = "2023-06-14T00:00:00"
        for field in core_fields:
            if field not in primary_header.keys():
                primary_header[field] = -99
    primary_header[BASE_NAME_KEY] = os.path.basename(path)
    primary_header[RAW_IMG_KEY] = path

    split_headers = tag_mef_extension_file_headers(
        primary_header=primary_header,
        extension_headers=split_headers,
        extension_key="BOARD_ID",
    )

    # Sometimes there are exptime keys
    for board_header in split_headers:
        if "EXPTIME" in board_header.keys():
            del board_header["EXPTIME"]

    return primary_header, split_data, split_headers


def load_winter_mef_image(
    path: str | Path,
) -> list[Image]:
    """
    Function to load winter mef images

    :param path: Path to image
    :return: list of images
    """
    images = open_mef_image(path, load_raw_winter_mef, extension_key="BOARD_ID")
    return images


def annotate_winter_subdet_headers(batch: ImageBatch) -> ImageBatch:
    """
    Annotate winter header with information on the subdetector

    :param batch: ImageBatch to annotate
    :return: ImageBatch where images have the updated header
    """
    new_batch = []
    for image in batch:
        data = image.get_data()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            _, med, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
            image["MEDCOUNT"] = med
            image["STDDEV"] = std

        subnx, subny, subnxtot, subnytot = (
            image["SUBNX"],
            image["SUBNY"],
            image["SUBNXTOT"],
            image["SUBNYTOT"],
        )

        mask = (
            (subdets["nx"] == subnx)
            & (subdets["ny"] == subny)
            & (subdets["nxtot"] == subnxtot)
            & (subdets["nytot"] == subnytot)
            & (subdets["boardid"] == image["BOARD_ID"])
        )
        assert np.sum(mask) == 1, (
            f"Subdet not found for nx={subnx}, ny={subny}, "
            f"nxtot={subnxtot}, nytot={subnytot} and boardid={image['BOARD_ID']}"
        )
        image["SUBDETID"] = int(subdets[mask]["subdetid"].iloc[0])
        image["RAWID"] = int(f"{image['EXPID']}_{str(image['SUBDETID']).rjust(2, '0')}")
        image["USTACKID"] = None

        if "DATASEC" in image.keys():
            del image["DATASEC"]

        # TODO: Write a little snippet to estimate the central RA/Dec from the pointing
        #  RA/Dec, BOARD_ID, SUBCOORD, and PA

        new_batch.append(image)
    new_batch = ImageBatch(new_batch)
    return new_batch


def get_raw_winter_mask(image: Image) -> np.ndarray:
    """
    Get mask for raw winter image.
    """
    data = image.get_data()
    header = image.header

    mask = np.zeros(data.shape)
    if header["BOARD_ID"] == 0:
        # Mask the outage in the bottom center
        mask[:500, 700:1600] = 1.0
        mask[1075:, :] = 1.0
        mask[:, 1950:] = 1.0
        mask[:20, :] = 1.0

    if header["BOARD_ID"] == 1:
        mask[:, 344:347] = 1.0
        mask[:, 998:1000] = 1.0
        mask[:, 1006:1008] = 1.0
        mask[260:262, :] = 1.0
        # Mask entire striped area to the right of the chip
        # mask[:, 1655:] = 1.0

        # Mask the low sensitivity regions around edges
        mask[1070:, :] = 1.0
        mask[:20, :] = 1.0
        mask[:, :75] = 1.0
        mask[:, 1961:] = 1.0

    if header["BOARD_ID"] == 2:
        mask[1060:, :] = 1.0
        mask[:, 1970:] = 1.0
        mask[:55, :] = 1.0
        mask[:, :20] = 1.0
        mask[:475, 406:419] = 1.0
        mask[350:352, :] = 1.0
        mask[260:287, :66] = 1.0
        mask[:, 1564:1567] = 1.0
        mask[:, 1931:] = 1.0

    if header["BOARD_ID"] == 3:
        mask[1085:, :] = 1.0
        mask[:, 1970:] = 1.0
        mask[:55, :] = 1.0
        mask[:, :20] = 1.0

        # Mask outages on top right and bottom right
        mask[:180, 1725:] = 1.0
        mask[1030:, 1800:] = 1.0

    if header["BOARD_ID"] == 4:
        # # Mask the region to the top left
        mask[610:, :250] = 1.0
        # # There seems to be a dead spot in the middle of the image
        mask[503:518, 390:405] = 1.0

        # Mask the edges with low sensitivity due to masking
        mask[:, 1948:] = 1.0
        mask[:, :61] = 1.0
        mask[:20, :] = 1.0
        mask[1060:, :] = 1.0

        # Mask a vertical strip
        mask[:, 998:1002] = 1.0

        # Mask the outage to the right
        mask[145:, 1735:] = 1.0
        # mask[data > 40000] = 1.0

        # Mask random vertical strip
        mask[:, 1080:1085] = 1.0

    if header["BOARD_ID"] == 5:
        # Mask the outage in the top-right.
        mask[700:, 1200:1900] = 1.0
        mask[1072:, :] = 1.0
        mask[:, 1940:] = 1.0
        mask[:15, :] = 1.0

    return mask.astype(bool)
