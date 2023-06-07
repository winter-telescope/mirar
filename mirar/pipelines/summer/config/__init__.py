"""
Module for summer-specific files
"""
import os

from mirar.pipelines.summer.config.constants import (
    DB_NAME,
    PIPELINE_NAME,
    SUMMER_GAIN,
    SUMMER_PIXEL_SCALE,
)
from mirar.pipelines.summer.config.schema import get_summer_schema_path
from mirar.processors.utils.cal_hunter import CalRequirement

summer_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(summer_dir, "files")
summer_mask_path = os.path.join(summer_dir, "mask.fits")
summer_weight_path = os.path.join(summer_dir, "weight.fits")

sextractor_astrometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "astrom.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "astrom.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}


sextractor_photometry_config = {
    "config_path": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "filter_path": os.path.join(astromatic_config_dir, "default.conv"),
    "parameter_path": os.path.join(astromatic_config_dir, "photom.param"),
    "starnnw_path": os.path.join(astromatic_config_dir, "default.nnw"),
}

sextractor_candidates_config = {
    "cand_det_sextractor_config": os.path.join(astromatic_config_dir, "photomCat.sex"),
    "cand_det_sextractor_nnw": os.path.join(astromatic_config_dir, "default.nnw"),
    "cand_det_sextractor_filter": os.path.join(astromatic_config_dir, "default.conv"),
    "cand_det_sextractor_params": os.path.join(astromatic_config_dir, "Scorr.param"),
}

scamp_path = os.path.join(astromatic_config_dir, "scamp.conf")

swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")

psfex_config_path = os.path.join(astromatic_config_dir, "photom.psfex")

summer_cal_requirements = [
    CalRequirement(
        target_name="bias", required_field="EXPTIME", required_values=["0.0"]
    ),
    CalRequirement(
        target_name="flat",
        required_field="FILTERID",
        required_values=["u", "g", "r", "i"],
    ),
]
