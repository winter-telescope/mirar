"""
Module for GMOS-specific constants
"""

# pylint: disable=duplicate-code

PIPELINE_NAME = "gmos"

GMOS_PIXEL_SCALE = 0.0807  # arcseconds per pixel
GMOS_ANET_INDEX_SERIES = "5200"

# This is what we support right now, not an exhaustive list
GMOS_FILTERS = ["u'", "g'", "r'", "i'", "z'", "y", "z"]
