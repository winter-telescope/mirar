from pathlib import Path

from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import OBSCLASS_KEY, TARGET_KEY, core_fields
from mirar.pipelines.spring.config.constants import SPRING_GAIN

# def load_raw_spring_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
#     """
#     Function to load a raw GIT image

#     :param path: path of file
#     :return: data and header of image


#     """
#     data, header = open_fits(path)

#     ## INSERT CODE HERE ##

#     for field in core_fields:
#         if field not in header:
#             raise KeyError(f"Core field {field} not found in header")


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
            header[OBSCLASS_KEY] = header["OBSTYPE"].strip().upper()
        else:
            header[OBSCLASS_KEY] = "SCIENCE"

    # -----------------------------
    # Target identification
    # -----------------------------
    if TARGET_KEY not in header:
        if "TARGNAME" in header and header["TARGNAME"] is not None:
            header[TARGET_KEY] = header["TARGNAME"].strip().upper()
        else:
            header[TARGET_KEY] = "UNKNOWN"

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
    header["RAWPATH"] = ""

    # -----------------------------
    # Final validation
    # -----------------------------
    for field in core_fields:
        if field not in header:
            raise KeyError(f"Core field {field} not found in header for {path}")

    return data, header


# def load_raw_spring_image(path: str | Path) -> Image:
#     """
#     Function to load a raw GIT image

#     :param path: Path to the raw image
#     :return: Image object
#     """
#     return open_raw_image(path, load_raw_spring_fits)


def load_raw_spring_image(path: str | Path) -> Image:
    return open_raw_image(path, load_raw_spring_fits)
