"""
Module for LMI-specific constants
"""

# pylint: disable=duplicate-code

PIPELINE_NAME = "lmi"
LMI_GAIN = 2.89  # in electron / ADU
LMI_READ_NOISE = 6  # in electron

LMI_PIXEL_SCALE = 0.24  # arcsec/pixel after 2x2 binning
LMI_NONLINEAR_LEVEL = 60000
LMI_FILTERS = ["u", "g", "r", "i", "z", "open"]

LMI_WIDTH_DEG = 12.3 / 60.0  # LMI is 12.3 arc minutes each side
