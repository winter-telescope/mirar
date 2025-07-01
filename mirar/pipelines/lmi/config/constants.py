"""
Module for LMI-specific constants
"""

# pylint: disable=duplicate-code

PIPELINE_NAME = "lmi"
LMI_GAIN = 2.89  # in electron / ADU
LMI_READ_NOISE = 6  # in electron

LMI_PIXEL_SCALE = 0.12  # arcsec/pixel
LMI_NONLINEAR_LEVEL = 100000
LMI_FILTERS = ["u", "g", "r", "i", "z", "open"]

LMI_WIDTH = 12.3 / 60.0  # LMI is 12.3 arc minutes each side
