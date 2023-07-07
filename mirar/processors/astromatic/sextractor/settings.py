"""
Module with sextractor utilities
"""
import os
from pathlib import Path

from mirar.processors.astromatic.config import astromatic_config_dir

default_param_path = os.path.join(astromatic_config_dir, "temp.param")
default_conv_path = os.path.join(astromatic_config_dir, "sex.conv")
default_config_path = os.path.join(astromatic_config_dir, "sex.config")
default_starnnw_path = os.path.join(astromatic_config_dir, "default.nnw")


def parse_sextractor_config(config_file: str | Path):
    """
    Parse a sextractor config file into a dictionary
    param config_file: path to sextractor config file
    """
    with open(config_file, "r") as f:
        data = f.readlines()

    keys, values = [], []
    config_dict = {}
    for row in data:
        if row == "\n":
            continue
        if row[0] == "#":
            continue
        key = row.split(" ")[0]
        if key == "":
            continue
        keys.append(key.split("\t")[0])
        for val in row.split(" ")[1:]:
            if val != "":
                if val == "#":
                    val = ""
                values.append(val)
                config_dict[key] = val.split("\t")[0]
                break
    return config_dict


def write_sextractor_config_to_file(config_dict: dict, config_filename: str | Path):
    """
    Write a sextractor config file from a dictionary
    param config_dict: dictionary of sextractor config parameters
    param config_filename: path to write config file
    """
    with open(config_filename, "w") as f:
        for key in config_dict:
            f.write(f"{key.ljust(20, ' ')}   {config_dict[key]} \n")


def write_param_file(param_path: str = default_param_path, params: list = None):
    """
    Write a default parameter file for sextractor
    param param_path: path to write parameter file
    param params: list of parameters to write. If None, will write default.
    """
    if params is None:
        params = [
            "X_IMAGE",
            "Y_IMAGE",
            "ALPHA_J2000",
            "DELTA_J2000",
            "MAG_AUTO",
            "MAGERR_AUTO",
            "ELLIPTICITY",
            "FWHM_IMAGE",
            "FLAGS",
        ]
    with open(param_path, "w", encoding="utf8") as param_f:
        for param in params:
            param_f.write(f"{param}\n")


def write_conv_file(conv_path: str = default_conv_path):
    """
    Write a default convolution file for sextractor
    """
    convol = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""
    with open(conv_path, "w", encoding="utf8") as conv_f:
        conv_f.write(convol)


def write_config_file(
    param_path: str = default_param_path,
    conv_path: str = default_conv_path,
    config_path: str = default_config_path,
    saturation_key: str = "SATURATE",
):
    """
    Write a default configuration file for sextractor
    """
    configs = f"""
#-------------------------------- Catalog ------------------------------------

CATALOG_NAME     temp.cat       # name of the output catalog
CATALOG_TYPE     ASCII_HEAD      # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
PARAMETERS_NAME  {param_path}     # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   5              # minimum number of pixels above threshold
DETECT_THRESH    3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

FILTER           Y              # apply filter for detection (Y or N)?
FILTER_NAME      {conv_path}       # name of the file containing the filter

DEBLEND_NTHRESH  16             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.02           # Minimum contrast parameter for deblending

CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency

MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                # NONE, BLANK or CORRECT

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                # <min_radius>



MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             0.0            # detector gain in e-/ADU
PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM      1.2            # stellar FWHM in arcsec
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_SIZE        64             # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  check.fits     # Filename for the check-image

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK  3000           # number of objects in stack
MEMORY_PIXSTACK  300000         # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE     QUIET          # can be QUIET, NORMAL or FULL
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output
SATUR_KEY        {saturation_key}       # keyword for saturation level (in ADUs)
"""
    with open(config_path, "w", encoding="utf8") as conf_f:
        conf_f.write(configs)
