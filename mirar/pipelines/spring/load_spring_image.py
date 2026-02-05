from pathlib import Path

from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import BASE_NAME_KEY, OBSCLASS_KEY, TARGET_KEY, core_fields
from mirar.pipelines.spring.config.constants import SPRING_GAIN


def load_raw_spring_fits(path: str | Path):
    data, header = open_fits(path)

    # -----------------------------
    # Time
    # -----------------------------
    # DATE-OBS already exists and is usable
    if "EXPTIME" not in header and "AEXPTIME" in header:
        header["EXPTIME"] = header["AEXPTIME"]

    # -----------------------------
    # Pointing (prefer degrees)
    # -----------------------------
    if "RA" not in header and "RADEG" in header:
        header["RA"] = header["RADEG"]

    if "DEC" not in header and "DECDEG" in header:
        header["DEC"] = header["DECDEG"]

    # -----------------------------
    # Instrument Identity
    # -----------------------------
    header.setdefault("INSTRUME", "SPRING")

    # -----------------------------
    # Filter handling
    # -----------------------------
    if "FILTER" not in header or not str(header["FILTER"]).strip():
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
    if "GAIN" not in header:
        header["GAIN"] = SPRING_GAIN

    # -----------------------------
    # Miscellaneous statement of header properties
    # -----------------------------
    header["COADDS"] = 1
    header["CALSTEPS"] = ""
    header["PROCFAIL"] = 1
    header["RAWPATH"] = path.as_posix()
    header[BASE_NAME_KEY] = Path(path).name

    # -----------------------------
    # Final validation
    # -----------------------------
    for field in core_fields:
        if field not in header:
            raise KeyError(f"Core field {field} not found in header for {path}")

    return data, header


def load_raw_spring_image(path: str | Path) -> Image:
    return open_raw_image(path, load_raw_spring_fits)
