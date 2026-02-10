from pathlib import Path

import numpy as np
from astropy.stats import sigma_clipped_stats

from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import BASE_NAME_KEY, OBSCLASS_KEY, TARGET_KEY, core_fields
from mirar.pipelines.spring.config.constants import SPRING_GAIN


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
    # -----------------------------
    # Final validation
    # -----------------------------
    for field in core_fields:
        if field not in header:
            raise KeyError(f"Core field {field} not found in header for {path}")

    data = data.astype("float32")
    return data, header


def load_raw_spring_image(path: str | Path) -> Image:
    return open_raw_image(path, load_raw_spring_fits)
