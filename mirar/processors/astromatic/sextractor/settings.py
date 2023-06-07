"""
Module with sextractor utilities
"""
import os

from mirar.processors.astromatic.config import astromatic_config_dir

default_param_path = os.path.join(astromatic_config_dir, "temp.param")
default_conv_path = os.path.join(astromatic_config_dir, "sex.conv")
default_config_path = os.path.join(astromatic_config_dir, "sex.config")
default_starnnw_path = os.path.join(astromatic_config_dir, "default.nnw")


def write_param_file(param_path: str = default_param_path):
    params = """X_IMAGE
Y_IMAGE
ALPHA_J2000
DELTA_J2000
MAG_AUTO
MAGERR_AUTO
ELLIPTICITY
FWHM_IMAGE
FLAGS"""
    with open(param_path, "w") as pf:
        pf.write(params)


def write_conv_file(conv_path: str = default_conv_path):
    convol = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""
    with open(conv_path, "w") as cf:
        cf.write(convol)


def write_config_file(
    param_path: str = default_param_path,
    conv_path: str = default_conv_path,
    config_path: str = default_config_path,
    saturation: float = 55000.0,
):
    configs = f"""
#-------------------------------- Catalog ------------------------------------

CATALOG_NAME     temp.cat       # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
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
SATUR_LEVEL      {saturation}        # level (in ADUs) at which arises saturation
"""
    with open(config_path, "w") as pf:
        pf.write(configs)
