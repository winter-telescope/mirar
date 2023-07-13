# Copyright 2015 Fred Moolekamp
# BSD 3-clause license
"""
Functions to convert FITS files or astropy Tables to FITS_LDAC files and
vice versa.
"""
import tempfile
import warnings
from pathlib import Path

import astropy.io
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning


def convert_hdu_to_ldac(
    hdu: astropy.io.fits.BinTableHDU | astropy.io.fits.TableHDU,
) -> tuple[astropy.io.fits.BinTableHDU, astropy.io.fits.BinTableHDU]:
    """
    Convert an hdu table to a fits_ldac table (format used by astromatic suite)

    Parameters
    ----------
    hdu: `astropy.io.fits.BinTableHDU` or `astropy.io.fits.TableHDU`
        HDUList to convert to fits_ldac HDUList

    Returns
    -------
    tbl1: `astropy.io.fits.BinTableHDU`
        Header info for fits table (LDAC_IMHEAD)
    tbl2: `astropy.io.fits.BinTableHDU`
        Data table (LDAC_OBJECTS)
    """
    tblhdr = np.array([hdu.header.tostring(",")])
    col1 = fits.Column(name="Field Header Card", array=tblhdr, format="13200A")
    cols = fits.ColDefs([col1])
    tbl1 = fits.BinTableHDU.from_columns(cols)
    tbl1.header["TDIM1"] = f"(80, {len(hdu.header)})"
    tbl1.header["EXTNAME"] = "LDAC_IMHEAD"
    tbl2 = fits.BinTableHDU(hdu.data)
    tbl2.header["EXTNAME"] = "LDAC_OBJECTS"
    return tbl1, tbl2


def convert_table_to_ldac(tbl: astropy.table.Table) -> astropy.io.fits.HDUList:
    """
    Convert an astropy table to a fits_ldac

    Parameters
    ----------
    tbl: `astropy.table.Table`
        Table to convert to ldac format
    Returns
    -------
    hdulist: `astropy.io.fits.HDUList`
        FITS_LDAC hdulist that can be read by astromatic software
    """
    table = tbl.copy()
    # Cannot save "object"-type fields via fits
    del_list = [x for x in table.dtype.names if table.dtype[x].kind == "O"]
    table.remove_columns(del_list)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        with tempfile.NamedTemporaryFile(suffix=".fits", mode="rb+") as temp_file:
            table.write(temp_file, format="fits")
            temp_file.seek(0)
            with fits.open(temp_file, mode="update") as hdulist:
                tbl1, tbl2 = convert_hdu_to_ldac(hdulist[1].copy())
                new_hdulist = [hdulist[0].copy(), tbl1, tbl2]
                new_hdulist = fits.HDUList(new_hdulist)
    return new_hdulist


def save_table_as_ldac(tbl: astropy.table.Table, file_path: str, **kwargs):
    """
    Save a table as a fits LDAC file

    Parameters
    ----------
    tbl: `astropy.table.Table`
        Table to save
    file_path: str
        Filename to save table
    kwargs:
        Keyword arguments to pass to hdulist.writeto
    """
    hdulist = convert_table_to_ldac(tbl)
    hdulist.writeto(file_path, **kwargs)


def get_table_from_ldac(file_path: str | Path, frame: int = 1) -> astropy.table.Table:
    """
    Load an astropy table from a fits_ldac by frame (Since the ldac format has column
    info for odd tables, giving it twce as many tables as a regular fits BinTableHDU,
    match the frame of a table to its corresponding frame in the ldac file).

    Parameters
    ----------
    file_path: str
        Name of the file to open
    frame: int
        Number of the frame in a regular fits file
    """
    if frame > 0:
        frame = frame * 2
    tbl = Table.read(file_path, hdu=frame)
    return tbl
