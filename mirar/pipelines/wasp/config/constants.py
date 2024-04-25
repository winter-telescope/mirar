"""
Module for WASP-specific constants
"""

# pylint: disable=duplicate-code

PIPELINE_NAME = "wasp"
WASP_GAIN = 5.9  # in electron / ADU
WASP_READ_NOISE = 5  # in electron

WASP_PIXEL_SCALE = 0.18  # arcsec/pixel
WASP_NONLINEAR_LEVEL = 50000
WASP_FILTERS = ["u'", "g'", "r'", "i'", "z'"]
