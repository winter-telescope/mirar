"""
Module for SEDMv2-specific constants
"""

PIPELINE_NAME = "sedmv2"
SEDMV2_GAIN = 1
SEDMV2_PIXEL_SCALE = 0.322  # arcsec/pixel # 0.3515625 newer?
# The 6 arcmin field needs quads below 4 arcmin, which only 5200 reaches here
SEDMV2_ANET_INDEX_SERIES = "5200"
DB_NAME = PIPELINE_NAME
