"""
Module containing SEDMv2-specific paths
"""

import os

from mirar.pipelines.git.config.constants import (
    DB_NAME,
    GIT_GAIN,
    GIT_PIXEL_SCALE,
    PIPELINE_NAME,
)

git_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(git_dir, "files")

sextractor_photometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "photom.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}

sextractor_PSF_photometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "photomPSF.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}

sextractor_astrometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "astrom.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "astrom.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}


sextractor_candidates_config = {
    "cand_det_sextractor_config": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "cand_det_sextractor_nnw": os.path.join(astromatic_config_dir, "default.nnw"),
    "cand_det_sextractor_filter": os.path.join(astromatic_config_dir, "default.conv"),
    "cand_det_sextractor_params": os.path.join(astromatic_config_dir, "Scorr.param"),
}

sextractor_reference_config = {
    "config_path": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "parameter_path": os.path.join(astromatic_config_dir, "photom.param"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}

scamp_path = os.path.join(astromatic_config_dir, "scamp.conf")

swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")

psfex_config_path = os.path.join(astromatic_config_dir, "photom.psfex")

psfex_sci_config_path = os.path.join(astromatic_config_dir, "photom_sci.psfex")
