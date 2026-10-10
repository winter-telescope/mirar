"""
Module for summer-specific files
"""

# pylint: disable=duplicate-code

from pathlib import Path

from mirar.pipelines.summer.config.constants import (
    DB_NAME,
    PIPELINE_NAME,
    SUMMER_GAIN,
    SUMMER_PIXEL_SCALE,
)
from mirar.processors.astromatic.sextractor.config_type import (
    SextractorCandidateConfig,
    SextractorConfig,
)
from mirar.processors.utils.cal_hunter import CalRequirement

summer_dir = Path(__file__).parent

astromatic_config_dir = summer_dir / "files"
summer_mask_path = summer_dir / "mask.fits"
summer_weight_path = summer_dir / "weight.fits"

sextractor_astrometry_config: SextractorConfig = {
    "config_path": astromatic_config_dir / "astrom.sex",
    "filter_path": astromatic_config_dir / "default.conv",
    "parameter_path": astromatic_config_dir / "astrom.param",
    "starnnw_path": astromatic_config_dir / "default.nnw",
}


sextractor_photometry_config: SextractorConfig = {
    "config_path": astromatic_config_dir / "photomCat.sex",
    "filter_path": astromatic_config_dir / "default.conv",
    "parameter_path": astromatic_config_dir / "photom.param",
    "starnnw_path": astromatic_config_dir / "default.nnw",
}

sextractor_candidates_config: SextractorCandidateConfig = {
    "cand_det_sextractor_config": astromatic_config_dir / "photomCat.sex",
    "cand_det_sextractor_nnw": astromatic_config_dir / "default.nnw",
    "cand_det_sextractor_filter": astromatic_config_dir / "default.conv",
    "cand_det_sextractor_params": astromatic_config_dir / "Scorr.param",
}

scamp_path = astromatic_config_dir / "scamp.conf"

swarp_config_path = astromatic_config_dir / "config.swarp"

psfex_config_path = astromatic_config_dir / "photom.psfex"

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
