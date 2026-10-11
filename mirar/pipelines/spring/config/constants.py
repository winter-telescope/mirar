"""
Module for SPRING-specific constants
"""

PIPELINE_NAME = "spring"
SPRING_GAIN = 4  # ~4-ish according to Nathan
SPRING_PIXEL_SCALE = 0.411  # arcsec/pixel -- obtained from nova.astrometry.net tests
# SPRING shares the WINTER machine, which holds the 4200 index files
SPRING_ANET_INDEX_SERIES = "4200"
DB_NAME = PIPELINE_NAME
