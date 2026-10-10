"""
Module containing SEDMv2-specific paths
"""

from pathlib import Path

from mirar.pipelines.sedmv2.config.constants import (
    DB_NAME,
    PIPELINE_NAME,
    SEDMV2_GAIN,
    SEDMV2_PIXEL_SCALE,
)
from mirar.processors.astromatic.sextractor.config_type import (
    SextractorCandidateConfig,
    SextractorConfig,
)
from mirar.processors.utils.cal_hunter import CalRequirement

sedmv2_dir = Path(__file__).parent

astromatic_config_dir = sedmv2_dir / "files"
sedmv2_mask_path = sedmv2_dir / "mask.fits"
sedmv2_weight_path = sedmv2_dir / "weight.fits"

sextractor_photometry_config: SextractorConfig = {
    "config_path": astromatic_config_dir / "photomCat.sex",
    "filter_path": astromatic_config_dir / "default.conv",
    "parameter_path": astromatic_config_dir / "photom.param",
    "starnnw_path": astromatic_config_dir / "default.nnw",
}

sextractor_PSF_photometry_config: SextractorConfig = {
    "config_path": astromatic_config_dir / "photomCat.sex",
    "filter_path": astromatic_config_dir / "default.conv",
    "parameter_path": astromatic_config_dir / "photomPSF.param",
    "starnnw_path": astromatic_config_dir / "default.nnw",
}

sextractor_astrometry_config: SextractorConfig = {
    "config_path": astromatic_config_dir / "astrom.sex",
    "filter_path": astromatic_config_dir / "default.conv",
    "parameter_path": astromatic_config_dir / "astrom.param",
    "starnnw_path": astromatic_config_dir / "default.nnw",
}


sextractor_candidate_config: SextractorCandidateConfig = {
    "cand_det_sextractor_config": astromatic_config_dir / "photomCat.sex",
    "cand_det_sextractor_nnw": astromatic_config_dir / "default.nnw",
    "cand_det_sextractor_filter": astromatic_config_dir / "default.conv",
    "cand_det_sextractor_params": astromatic_config_dir / "Scorr.param",
}

sextractor_reference_config: SextractorConfig = {
    "config_path": astromatic_config_dir / "photomCat.sex",
    "parameter_path": astromatic_config_dir / "photom.param",
    "filter_path": astromatic_config_dir / "default.conv",
    "starnnw_path": astromatic_config_dir / "default.nnw",
}

scamp_path = astromatic_config_dir / "scamp.conf"

swarp_config_path = astromatic_config_dir / "config.swarp"

psfex_config_path = astromatic_config_dir / "photom.psfex"

sedmv2_cal_requirements = [
    CalRequirement(
        target_name="flat",
        required_field="FILTERID",
        required_values=["u", "g", "r", "i"],
    ),
]
