"""
Module for LRIS-specific constants
"""

# pylint: disable=duplicate-code

# https://www2.keck.hawaii.edu/observing/kecktelgde/ktelinstupdate.pdf
PIPELINE_NAME = "lris"
LRIS_GAIN = 1.6  # in electron / ADU
LRIS_READ_NOISE = 2.5  # in electron

LRIS_PIXEL_SCALE = 0.135  # arcsec/pixel
LRIS_NONLINEAR_LEVEL = 100000
LRIS_FILTERS = ["u'"]

LRIS_WIDTH = 7.8 / 60.0  # LRIS is 6 x 7.8 arc minutes
