import logging
import warnings
from pathlib import Path

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import BASE_NAME_KEY, OBSCLASS_KEY, TARGET_KEY
from mirar.pipelines.spring.config.constants import SPRING_GAIN
from mirar.pipelines.spring.constants import imgtype_dict, spring_filters_map
from mirar.pipelines.spring.models import default_program, itid_dict

logger = logging.getLogger(__name__)


def load_raw_spring_fits(path: str | Path):
    data, header = open_fits(path)

    header.remove("BITPIX", ignore_missing=True)
    header.remove("BZERO", ignore_missing=True)
    header.remove("BSCALE", ignore_missing=True)
    # -----------------------------
    # Time
    # -----------------------------
    # DATE-OBS already exists and is usable
    if "EXPTIME" not in header and "AEXPTIME" in header:
        header["EXPTIME"] = header["AEXPTIME"]

    # -----------------------------
    # Pointing (prefer degrees)
    # -----------------------------
    header["RA"] = header["RADEG"]
    header["DEC"] = header["DECDEG"]

    # -----------------------------
    # Instrument Identity
    # -----------------------------
    header.setdefault("INSTRUME", "SPRING")

    # -----------------------------
    # Filter handling
    # -----------------------------
    if "FILTER" not in header or not str(header["FILTER"]).strip():
        header["EPSPP"] = (
            sigma_clipped_stats(data, sigma=5)[1] * SPRING_GAIN / header["EXPTIME"]
        )
        if header["EPSPP"] <= 100:
            header["FILTER"] = "Y"
        elif header["EPSPP"] >= 200 and header["EPSPP"] <= 800:
            header["FILTER"] = "J"
        elif header["EPSPP"] >= 1500 and header["EPSPP"] <= 4000:
            header["FILTER"] = "H"
        elif header["EPSPP"] >= 8000:
            header["FILTER"] = "K"
        else:
            header["FILTER"] = "UNKNOWN"

    if header["FILTER"] in spring_filters_map:
        header["FID"] = int(spring_filters_map[header["FILTER"]])
    else:
        header["FID"] = -99

    # -----------------------------
    # Observation classification
    # -----------------------------
    if OBSCLASS_KEY not in header:
        if "OBSTYPE" in header and header["OBSTYPE"] is not None:
            header[OBSCLASS_KEY] = header["OBSTYPE"].strip().lower()
        else:
            header[OBSCLASS_KEY] = "science"

    # -----------------------------
    # Target identification, same logic as WINTER.
    # -----------------------------
    target = f"field_{header['FIELDID']}"
    if ("SCHDNAME" in header.keys()) & ("OBHISTID" in header.keys()):
        if header["SCHDNAME"] != "":
            target = f"{header['SCHDNAME']}_{header['OBHISTID']}"
    elif TARGET_KEY in header.keys():
        if header[TARGET_KEY] != "":
            target = header[TARGET_KEY]

    if header[OBSCLASS_KEY].lower() in [
        "dark",
        "bias",
        "focus",
        "pointing",
        "flat",
        "test",
        "corrupted",
    ]:
        target = header[OBSCLASS_KEY].lower()
    header[TARGET_KEY] = target
    # -----------------------------
    # Camera GAIN
    # -----------------------------
    # if "GAIN" not in header:
    header["GAIN"] = SPRING_GAIN

    # -----------------------------
    # Miscellaneous statement of header properties
    # -----------------------------
    header["COADDS"] = 1
    header["CALSTEPS"] = ""
    header["PROCFAIL"] = 1
    header["RAWPATH"] = path.as_posix()
    header[BASE_NAME_KEY] = Path(path).name
    header["MEDCOUNT"] = np.nanmedian(data)

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

    header["UTCTIME"] = Time(header["UTCISO"], format="iso").isot
    date_t = Time(header["UTCTIME"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        header["NIGHTDATE"] = date_t.to_datetime().strftime("%Y-%m-%d")

    header["IMGTYPE"] = header[OBSCLASS_KEY]
    if header["IMGTYPE"] == "test":
        header["IMGTYPE"] = "other"

    if not header["IMGTYPE"] in imgtype_dict:
        logger.error(f"Unknown image type: {header['IMGTYPE']}")
        header["ITID"] = itid_dict[imgtype_dict["corrupted"]]
    else:
        header["ITID"] = itid_dict[imgtype_dict[header["IMGTYPE"]]]

    header["EXPMJD"] = header["MJD-OBS"]

    header["RAWID"] = int((date_t.mjd - 59000.0) * 86400.0)  # seconds since 59000 MJD
    data = data.astype("float32")
    return data, header


def load_raw_spring_image(path: str | Path) -> Image:
    return open_raw_image(path, load_raw_spring_fits)
