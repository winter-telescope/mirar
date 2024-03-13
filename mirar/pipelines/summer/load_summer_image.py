"""
Module with functions to load raw and processed summer images
"""

import logging
import warnings
from pathlib import Path

import astropy
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import (
    BASE_NAME_KEY,
    GAIN_KEY,
    LATEST_SAVE_KEY,
    OBSCLASS_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    TARGET_KEY,
    __version__,
)
from mirar.pipelines.summer.models import DEFAULT_FIELD, SUMMER_NIGHT_FORMAT

logger = logging.getLogger(__name__)


def load_raw_summer_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw summer image and add/modify the required headers
    Args:
        path: Path to the raw image

    Returns: [image data, image header]

    """
    if isinstance(path, str):
        path = Path(path)
    data, header = open_fits(path)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        header[OBSCLASS_KEY] = header["OBSTYPE"].lower()
        header["UTCTIME"] = str(header["UTCSHUT"]).replace(" ", "T")

        header["MJD-OBS"] = header["OBSMJD"]
        header["DATE-OBS"] = Time(header["OBSMJD"], format="mjd").isot

        # If the image is a calibration image, set the target name to the OBSTYPE.
        # If it is a science image, it is either a field observation, or a ToO with
        # a target name. If it is a field observation, set the target name to the
        # field ID. If it is a ToO, set the target name to the TARGNAME.
        if header[OBSCLASS_KEY] == "science":
            if "TARGNAME" in header:
                target_name = header["TARGNAME"]
            else:
                target_name = header[OBSCLASS_KEY]
            if target_name == "":
                target_name = f"field_{header['FIELDID']}"
        else:
            target_name = header[OBSCLASS_KEY]

        try:
            header[TARGET_KEY] = target_name.lower()
        except (ValueError, AttributeError):
            header[TARGET_KEY] = target_name

        crd = SkyCoord(ra=header["RA"], dec=header["DEC"], unit=(u.deg, u.deg))
        header["RA"] = crd.ra.deg
        header["DEC"] = crd.dec.deg

        header["CRVAL1"] = header["RA"]
        header["CRVAL2"] = header["DEC"]

        tel_crd = SkyCoord(
            ra=header["TELRA"],
            dec=header["TELDEC"],
            unit=(u.deg, u.deg),
        )
        header["TELRA"] = tel_crd.ra.deg
        header["TELDEC"] = tel_crd.dec.deg
        header["BZERO"] = 0

        header[LATEST_SAVE_KEY] = path.as_posix()
        header[RAW_IMG_KEY] = path.as_posix()

        data = data * 1.0  # pylint: disable=no-member

        if "other" in header["FILTERID"]:
            header["FILTERID"] = "r"

        header[PROC_HISTORY_KEY] = ""

        base_name = path.name
        header[BASE_NAME_KEY] = base_name

        pipeline_version = __version__
        pipeline_version_padded_str = "".join(
            [x.rjust(2, "0") for x in pipeline_version.split(".")]
        )

        obstime = Time(header["UTCISO"], format="iso")
        header["EXPID"] = int(
            (obstime.mjd - 59000.0) * 86400.0
        )  # seconds since 60000 MJD

        header["PROCID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))

        header["OBSDATE"] = int(header["UTC"].split("_")[0])

        header["TIMEUTC"] = header["UTCISO"]

        header["NIGHTDATE"] = obstime.to_datetime().strftime(SUMMER_NIGHT_FORMAT)
        header["EXPMJD"] = header["OBSMJD"]
        header["OBS-DATE"] = Time(header["OBSMJD"], format="mjd").isot

        default_id = 0

        for key in ["PROGID", "OBSID"]:
            if key not in header.keys():
                # logger.warning(f"No {key} found in header of {path}")
                header[key] = default_id
            else:
                try:
                    header[key] = int(header[key])
                except ValueError:
                    header[key] = default_id

        if "SUBPROG" not in header.keys():
            # logger.warning(f"No SUBPROG found in header of {path}")
            header["SUBPROG"] = "none"

        header["FILTER"] = header["FILTERID"]
        header["FID"] = header["FILPOS"]
        try:
            header["SHUTOPEN"] = Time(header["SHUTOPEN"], format="iso").jd
        except ValueError:
            logger.debug(f"Error parsing 'SHUTOPEN' of {path}: ({header['SHUTOPEN']})")

        try:
            header["SHUTCLSD"] = Time(header["SHUTCLSD"], format="iso").jd
        except ValueError:
            logger.debug(f"Error parsing 'SHUTCLSD' of {path}: ({header['SHUTCLSD']})")

        header["PROCFLAG"] = 0

        header[PROC_FAIL_KEY] = ""

        sunmoon_keywords = [
            "MOONRA",
            "MOONDEC",
            "MOONILLF",
            "MOONPHAS",
            "MOONALT",
            "SUNAZ",
            "SUNALT",
        ]
        for key in sunmoon_keywords:
            val = 0
            if key in header.keys():
                if header[key] not in [""]:
                    val = header[key]
            header[key] = val

        itid_dict = {
            "SCIENCE": 1,
            "BIAS": 2,
            "FLAT": 2,
            "DARK": 2,
            "FOCUS": 3,
            "POINTING": 4,
            "OTHER": 5,
        }

        if not header[OBSCLASS_KEY] in itid_dict:
            header["ITID"] = 5
        else:
            header["ITID"] = itid_dict[header[OBSCLASS_KEY]]

        if header["FIELDID"] == "radec":
            header["FIELDID"] = DEFAULT_FIELD

        if header["ITID"] != 1:
            header["FIELDID"] = DEFAULT_FIELD

        if "COADDS" not in header.keys():
            header["COADDS"] = 1

        if "PROGPI" not in header.keys():
            header["PROGPI"] = "?"

        if header["PROGID"] in ["", "WINTER"]:
            header["PROGID"] = 2

        try:
            header["PROGID"] = int(header["PROGID"])
        except ValueError:
            try:
                progpi = header["PROGID"]
                header["PROGID"] = int(header["PROGPI"])
                header["PROGPI"] = progpi
            except KeyError:
                header["PROGID"] = 1

        # TODO Figure out how if database query is required for this.
        header["PUID"] = header["PROGID"]
        crds = SkyCoord(ra=header["RA"], dec=header["DEC"], unit=(u.deg, u.deg))
        header["RA"] = crds.ra.deg
        header["DEC"] = crds.dec.deg

        # TODO Write a processor to calculate the split QID
        header["QID"] = 1
        header["PROCSTATUS"] = 0

        # TODO Figure out what to do about primary keys
        header["RAWID"] = header["EXPID"]  # + subdetid

        if GAIN_KEY not in header.keys():
            header[GAIN_KEY] = 1.0

    return data, header  # pylint: disable=no-member


def load_raw_summer_image(path: str | Path) -> Image:
    """
    Function to load a raw summer image and add/modify the required headers

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_summer_fits)


def load_proc_summer_image(path: str) -> Image:
    """
    Function to load a processed summer image and add/modify the required headers
    Args:
        path: Path to the processed image

    Returns: [image data, image header]

    """
    data, header = open_fits(path)

    if "ZP" not in header.keys():
        header["ZP"] = header["ZP_AUTO"]
        header["ZP_std"] = header["ZP_AUTO_std"]

    header["CENTRA"] = header["CRVAL1"]
    header["CENTDEC"] = header["CRVAL2"]

    if "TARGET" in header.keys():
        header[TARGET_KEY] = header["TARGET"]

    pipeline_version = __version__
    pipeline_version_padded_str = "".join(
        [x.rjust(2, "0") for x in pipeline_version.split(".")]
    )
    header["DIFFID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))

    if "PROCID" not in header.keys():
        header["PROCID"] = header["DIFFID"]

    header["TIMEUTC"] = header["UTCISO"]
    data[data == 0] = np.nan
    return Image(data, header)
