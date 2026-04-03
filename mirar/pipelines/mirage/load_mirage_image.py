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
from mirar.pipelines.mirage.config.constants import MIRAGE_GAIN
from mirar.pipelines.mirage.constants import imgtype_dict, mirage_filters_map
from mirar.pipelines.mirage.models import default_program, itid_dict

logger = logging.getLogger(__name__)


def load_raw_mirage_fits(path: str | Path):
    data, header = open_fits(path)

    header.remove("BITPIX", ignore_missing=True)
    header.remove("BZERO", ignore_missing=True)
    header.remove("BSCALE", ignore_missing=True)
    # -----------------------------
    # Time
    # -----------------------------
    # DATE-OBS already exists and is usable
    # -----------------------------
    # Pointing (prefer degrees)
    # -----------------------------
    header["RA"] = header["CRVAL1"]
    header["DEC"] = header["CRVAL2"]
    header["RADEG"] = header["CRVAL1"]
    header["DECDEG"] = header["CRVAL2"]

    # -----------------------------
    # Instrument Identity
    # -----------------------------
    header.setdefault("INSTRUME", "MIRAGE")

    # -----------------------------
    # Filter handling
    # -----------------------------
    if "FILTER" not in header:
        object = header.get("OBJECT", "").strip()
        if "_" not in object:
            header["FILTER"] = "J"
        else:
            filter_name = object.split("_")[-1]
            header["FILTER"] = filter_name

    # -----------------------------
    # Observation classification
    # -----------------------------
    if OBSCLASS_KEY not in header:
        if header["OBJECT"].strip().lower() in ["scicam"]:
            header[OBSCLASS_KEY] = "flat"
        else:
            header[OBSCLASS_KEY] = "science"

    # -----------------------------
    # Target identification, same logic as WINTER.
    # -----------------------------
    header[TARGET_KEY] = header["OBJECT"].strip()

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
    header["GAIN"] = MIRAGE_GAIN

    # -----------------------------
    # Miscellaneous statement of header properties
    # -----------------------------
    header["COADDS"] = 1
    header["CALSTEPS"] = ""
    header["PROCFAIL"] = 1
    header["RAWPATH"] = path.as_posix()
    header[BASE_NAME_KEY] = Path(path).name
    header["MEDCOUNT"] = np.nanmedian(data)

    date_str = path.as_posix().split('/')[-1].split(".fits")[0].split("scicam_")[
        -1
    ]  # eg 20260401T123456
    date_t = Time(
        f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}T{date_str[9:11]}:{date_str[11:13]}:{date_str[13:15]}"
    )

    header['DATE-OBS'] = date_t.isot
    header["EXPMJD"] = date_t.mjd
    header["SAVEPATH"] = path.as_posix()
    data = data.astype("float32")
    return data, header


def load_raw_mirage_image(path: str | Path) -> Image:
    return open_raw_image(path, load_raw_mirage_fits)


def load_mirage_stack(
    path: str | Path,
) -> Image:
    """
    Load proc image


    :param path: Path to image
    :return: data and header
    """

    logger.debug(f"Loading {path}")
    data, header = open_fits(path)

    dirname = path.split("/mirage/")[0] + "/mirage/"
    wghtpath = header["WGHTPATH"]
    weight_pathname = wghtpath.split("/mirage/")[-1]
    new_weightpath = Path(dirname) / weight_pathname
    header["WGHTPATH"] = new_weightpath.as_posix()
    header["SAVEPATH"] = path

    if "PSFCAT" in header.keys():
        new_psfpath = Path(dirname) / header["PSFCAT"].split("/mirage/")[-1]
        header["PSFCAT"] = new_psfpath.as_posix()

    if "RFCTPATH" in header.keys():
        new_catpath = Path(dirname) / header["RFCTPATH"].split("/winter/")[-1]
        header["RFCTPATH"] = new_catpath.as_posix()

    if TARGET_KEY not in header.keys():
        if "TARGNAME" in header.keys():
            header[TARGET_KEY] = header["TARGNAME"]
    return Image(data=data, header=header)
