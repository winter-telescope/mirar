"""
Module containing WIRC-specific paths
"""

from pathlib import Path

from fastavro.schema import load_schema

from mirar.processors.astromatic.sextractor.config_type import (
    SextractorCandidateConfig,
    SextractorConfig,
)

wirc_file_dir = Path(__file__).parent / "files"
wirc_schema_dir = Path(__file__).parent / "schema"

wirc_mask_path = wirc_file_dir / "wirc_bad_feb_2018.fits"

sextractor_astrometry_config: SextractorConfig = {
    "config_path": wirc_file_dir / "stack.sex",
    "filter_path": wirc_file_dir / "default.conv",
    "parameter_path": wirc_file_dir / "astrom.param",
    "starnnw_path": wirc_file_dir / "default.nnw",
}

sextractor_photometry_config: SextractorConfig = {
    "config_path": wirc_file_dir / "matchcat.sex",
    "filter_path": wirc_file_dir / "default.conv",
    "parameter_path": wirc_file_dir / "astrom.param",
    "starnnw_path": wirc_file_dir / "default.nnw",
}

scamp_fp_path = wirc_file_dir / "scamp_fp.conf"

swarp_sp_path = wirc_file_dir / "second_pass.swarp"

sextractor_reference_config: SextractorConfig = {
    "config_path": wirc_file_dir / "photomCat.sex",
    "parameter_path": wirc_file_dir / "photom.param",
    "filter_path": wirc_file_dir / "default.conv",
    "starnnw_path": wirc_file_dir / "default.nnw",
}

sextractor_candidate_config: SextractorCandidateConfig = {
    "cand_det_sextractor_config": wirc_file_dir / "photomCat.sex",
    "cand_det_sextractor_nnw": wirc_file_dir / "default.nnw",
    "cand_det_sextractor_filter": wirc_file_dir / "default.conv",
    "cand_det_sextractor_params": wirc_file_dir / "Scorr.param",
}

psfex_path = wirc_file_dir / "photom.psfex"
