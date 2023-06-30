"""
Module component of autoastrometry, handling I/O
"""
import logging
import os
from pathlib import Path
from typing import Optional

import ephem
import numpy as np
from astropy.io import fits

from mirar.processors.astrometry.autoastrometry.sources import (
    BaseSource,
    SextractorSource,
)
from mirar.processors.astrometry.autoastrometry.utils import dec_str_2_deg, ra_str_2_deg

logger = logging.getLogger(__name__)


def parse_header(
    file_path: str | Path,
    temp_path: str | Path,
    pixel_scale: Optional[float] = None,
    pa: Optional[float] = None,
    inv: bool = False,
    user_ra_deg: Optional[float] = None,
    user_dec_deg: Optional[float] = None,
) -> tuple[int, int, float, float, float, float, float, float, float, float]:
    """
    Extract relevant info for astrometry from header of file

    :param file_path: File to parse
    :param temp_path: Temp path to write updated version to
    :param pixel_scale: Pixel scale
    :param pa: PA
    :param inv:
    :param user_ra_deg:
    :param user_dec_deg:
    :return: nxpix, nypix, cd11, cd12, cd21, cd22, crpix1, crpix2, cra, cdec
    """
    sci_ext = 0

    # Get some basic info from the header
    with fits.open(file_path) as hdu:
        hdu.verify("silentfix")

        header = hdu[sci_ext].header

        if np.logical_and(pixel_scale is not None, pa is None):
            pa = 0

        # Check for old-style WCS header
        if pixel_scale is None:
            old_wcs_type = False

            for hkey in header.keys():
                if hkey in ["CDELT1", "CDELT2"]:
                    old_wcs_type = True

            if old_wcs_type:
                key = "CDELT1"
                cdelt1 = header[key]
                key = "CDELT2"
                cdelt2 = header[key]

                try:
                    c_rot = 0
                    key = "CROTA1"
                    c_rot = header[key]
                    key = "CROTA2"
                    c_rot = header[key]
                except KeyError:
                    pass

                if (
                    np.sqrt(cdelt1**2 + cdelt2**2) < 0.1
                ):  # some images use CDELT to indicate nonstandard things
                    header["CD1_1"] = cdelt1 * np.cos(c_rot * np.pi / 180.0)
                    header["CD1_2"] = -cdelt2 * np.sin(c_rot * np.pi / 180.0)
                    header["CD2_1"] = cdelt1 * np.sin(c_rot * np.pi / 180.0)
                    header["CD2_2"] = cdelt2 * np.cos(c_rot * np.pi / 180.0)

        if np.logical_and(pixel_scale is not None, pa is not None):
            # Create WCS header information if pixel scale is specified
            pa_rad = pa * np.pi / 180.0
            px_scale_deg = pixel_scale / 3600.0

            if inv > 0:
                parity = -1
            else:
                parity = 1

            if user_ra_deg is not None:
                ra = user_ra_deg
            else:
                ra = ra_str_2_deg(header["CRVAL1"])

            if user_dec_deg is not None:
                dec = user_dec_deg
            else:
                dec = dec_str_2_deg(header["CRVAL2"])

            try:
                epoch = float(header.get("EPOCH", 2000))
            except KeyError:
                logger.warning("No EPOCH found in header. Assuming 2000")
                epoch = 2000.0

            try:
                equinox = float(
                    header.get("EQUINOX", epoch)
                )  # If RA and DEC are not J2000 then convert
            except KeyError:
                logger.warning("No EQUINOX found in header. Assuming 2000")
                equinox = 2000.0  # could be 'J2000'; try to strip off first character?

            if abs(equinox - 2000) > 0.5:
                logger.debug(f"Converting equinox from {equinox} to J2000")
                j2000 = ephem.Equatorial(
                    ephem.Equatorial(str(ra / 15), str(dec), epoch=str(equinox)),
                    epoch=ephem.J2000,
                )
                [ra, dec] = [ra_str_2_deg(j2000.ra), dec_str_2_deg(j2000.dec)]

            header["CD1_1"] = px_scale_deg * np.cos(pa_rad) * parity
            header["CD1_2"] = px_scale_deg * np.sin(pa_rad)
            header["CD2_1"] = -px_scale_deg * np.sin(pa_rad) * parity
            header["CD2_2"] = px_scale_deg * np.cos(pa_rad)
            header["CRPIX1"] = header["NAXIS1"] / 2
            header["CRPIX2"] = header["NAXIS2"] / 2
            header["CRVAL1"] = ra
            header["CRVAL2"] = dec
            header["CTYPE1"] = "RA---TAN"
            header["CTYPE2"] = "DEC--TAN"
            header["EQUINOX"] = 2000.0

            hdu[sci_ext].header = header
            hdu.writeto(
                temp_path, output_verify="silentfix", overwrite=True
            )  # ,clobber=True

        # Read the header info from the file.
        try:
            # no longer drawing RA and DEC from here.
            key = "NAXIS1"
            nxpix = header[key]
            key = "NAXIS2"
            nypix = header[key]
        except KeyError:
            err = f"Cannot find necessary WCS header keyword {key}"
            logger.debug(err)
            raise

        try:
            key = "CRVAL1"
            cra = float(header[key])
            key = "CRVAL2"
            cdec = float(header[key])

            key = "CRPIX1"
            crpix1 = float(header[key])
            key = "CRPIX2"
            crpix2 = float(header[key])

            key = "CD1_1"
            cd11 = float(header[key])
            key = "CD2_2"
            cd22 = float(header[key])
            key = "CD1_2"
            cd12 = float(header[key])  # deg / pix
            key = "CD2_1"
            cd21 = float(header[key])

            equinox = float(header.get("EQUINOX", 2000.0))
            if abs(equinox - 2000.0) > 0.2:
                logger.debug("Warning: EQUINOX is not 2000.0")

        except KeyError:
            if pixel_scale == -1:
                err = (
                    f"Cannot find necessary WCS header keyword '{key}' \n "
                    f"Must specify pixel scale (-px VAL) or "
                    f"provide provisional basic WCS info via CD matrix."
                )
                logger.error(err)
                raise
                # Some images might use CROT parameters, could try to be compatible with
                # this too...?

        # Wipe nonstandard hdu info from the header (otherwise this will confuse
        # verification)
        header_keys = list(header.keys())
        ctype_change = 0
        iraf_keys = []
        high_keys = []
        old_keys = []
        distortion_keys = []
        for hkey in header_keys:
            if (
                hkey == "RADECSYS"
                or hkey == "WCSDIM"
                or hkey.find("WAT") == 0
                or hkey.find("LTV") >= 0
                or hkey.find("LTM") == 0
            ):
                del header[hkey]
                iraf_keys.append(hkey)

            if (
                hkey.find("CO1_") == 0
                or hkey.find("CO2_") == 0
                or hkey.find("PV1_") == 0
                or hkey.find("PV2_") == 0
                or hkey.find("PC00") == 0
            ):
                del header[hkey]
                high_keys.append(hkey)

            if (
                hkey.find("CDELT1") == 0
                or hkey.find("CDELT2") == 0
                or hkey.find("CROTA1") == 0
                or hkey.find("CROTA2") == 0
            ):
                del header[hkey]
                old_keys.append(hkey)

            if (
                hkey.find("A_") == 0
                or hkey.find("B_") == 0
                or hkey.find("AP_") == 0
                or hkey.find("BP_") == 0
            ):
                del header[hkey]
                distortion_keys.append(hkey)

        if header["CTYPE1"] != "RA---TAN":
            logger.debug(f"Changing CTYPE1 from '{header['CTYPE1']}' to 'RA---TAN'")
            header["CTYPE1"] = "RA---TAN"
            ctype_change = 1

        if header["CTYPE2"] != "DEC--TAN":
            if ctype_change:
                logger.debug(f"Changing CTYPE2 from '{header['CTYPE2']}' to 'DEC--TAN'")
            header["CTYPE2"] = "DEC--TAN"
            ctype_change = 1

        wcs_key_check = [
            "CRVAL1",
            "CRVAL2",
            "CRPIX1",
            "CRPIX2",
            "CD1_1",
            "CD1_2",
            "CD2_2",
            "CD2_1",
            "EQUINOX",
            "EPOCH",
        ]

        header_format_change = False

        for w_key in wcs_key_check:
            if isinstance(w_key, str):
                try:
                    header[w_key] = float(header[w_key])
                    header_format_change = True
                except KeyError:
                    pass

        if len(iraf_keys) > 0:
            logger.warning(f"Removed nonstandard WCS keywords: {iraf_keys}")
        if len(high_keys) > 0:
            logger.warning(f"Removed higher-order WCS keywords: {high_keys}")
        if len(old_keys) > 0:
            logger.warning(f"Removed old-style WCS keywords: {old_keys}")
        if len(distortion_keys) > 0:
            logger.warning(f"Removed distortion WCS keywords: {distortion_keys}")

        if (
            len(high_keys) + len(distortion_keys) + ctype_change + header_format_change
            > 0
        ):
            # Rewrite and reload the image if the header
            # was modified in a significant way so
            # sextractor sees the same thing that we do.
            hdu[sci_ext].header = header
            hdu.writeto(
                temp_path, output_verify="silentfix", overwrite=True
            )  # ,clobber=True

    return nxpix, nypix, cd11, cd12, cd21, cd22, crpix1, crpix2, cra, cdec


def write_text_file(file_path: str, src_list: list[BaseSource]):
    """
    Write a text file with a list of sources

    :param file_path: Output file
    :param src_list: List of sources
    :return: None
    """
    logger.debug(f"Saving text file to {file_path}")

    with open(file_path, "w", encoding="utf8") as out:
        for src in src_list:
            out.write(f"{src.ra_deg:11.7f} {src.dec_deg:11.7f} {src.mag:5.2f}\n")


def write_region_file(
    file_path: str,
    src_list: list[BaseSource],
    color: str = "green",
    system: Optional[str] = None,
):
    """
    Write a region file

    :param file_path: Output path
    :param src_list: List of sources
    :param color: Colour to use
    :param system: system to use (default wcs)
    :return:
    """
    if system is None:
        system = "wcs"

    if system not in ["wcs", "img"]:
        err = f"Did not recognise system '{system}'. Valid values are 'wcs' and 'img'."
        logger.error(err)
        raise ValueError(err)

    logger.debug(f"Saving region file to {file_path}")

    with open(file_path, "w", encoding="utf8") as out:
        out.write(
            f"# Region file format: DS9 version 4.0\n"
            f'global color={color} font="helvetica 10 normal" select=1 '
            f"highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n"
        )

        if system == "wcs":
            out.write("fk5\n")
            for i, src in enumerate(src_list):
                out.write(
                    f"point({src.ra_deg:.7f},{src.dec_deg:.7f}) "
                    f"# point=boxcircle text={{{i + 1}}}\n"
                )
        elif system == "img":
            out.write("image\n")
            for i, src in enumerate(src_list):
                out.write(
                    f"point({src.x:.3f},{src.y:.3f}) "
                    f"# point=boxcircle text={{{i + 1}}}\n"
                )


def export_src_lists(
    img_src_list: list[SextractorSource],
    ref_src_list: list[BaseSource],
    base_output_path: str,
):
    """
    Export both image and reference source lists to txt/region files

    :param img_src_list: Image sources
    :param ref_src_list: Reference sources
    :param base_output_path: base output path
    :return: None
    """
    write_text_file(
        file_path=os.path.splitext(base_output_path)[0] + ".det.init.txt",
        src_list=img_src_list,
    )
    write_region_file(
        file_path=os.path.splitext(base_output_path)[0] + ".det.im.reg",
        src_list=img_src_list,
        color="red",
        system="img",
    )
    write_text_file(
        file_path=os.path.splitext(base_output_path)[0] + ".cat.txt",
        src_list=ref_src_list,
    )
    write_region_file(
        file_path=os.path.splitext(base_output_path)[0] + ".cat.wcs.reg",
        src_list=ref_src_list,
        color="green",
        system="wcs",
    )
