"""
Module containing GMOS-specific paths
"""

# pylint: disable=duplicate-code

import os

from mirar.pipelines.gmos.config.constants import PIPELINE_NAME

gmos_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(gmos_dir, "files")

sextractor_photometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "photom.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}

sextractor_astrometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "astrom.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "astrom.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}

scamp_path = os.path.join(astromatic_config_dir, "scamp.conf")

swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")

psfex_sci_config_path = os.path.join(astromatic_config_dir, "photom_sci.psfex")
