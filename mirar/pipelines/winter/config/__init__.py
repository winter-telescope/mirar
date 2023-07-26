"""
This file contains the configuration for the winter pipeline.
"""
from pathlib import Path

from mirar.processors.utils.cal_hunter import CalRequirement

PIPELINE_NAME = "winter"

winter_file_dir = Path(__file__).parent.joinpath("files")

sextractor_astrometry_config = {
    "config_path": winter_file_dir.joinpath("astrom.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("astrom.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}

sextractor_astromstats_config = sextractor_astrometry_config
sextractor_astromstats_config["parameter_path"] = winter_file_dir.joinpath(
    "astromstats.param"
)

sextractor_photometry_config = {
    "config_path": winter_file_dir.joinpath("photomCat.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("photom.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}


sextractor_anet_config = {
    "config_path": winter_file_dir.joinpath("astrom_anet.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("astrom.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}

sextractor_reference_config = {
    "config_path": winter_file_dir.joinpath("photomCat.sex"),
    "parameter_path": winter_file_dir.joinpath("photom.param"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}

sextractor_candidate_config = {
    "cand_det_sextractor_config": winter_file_dir.joinpath("photomCat.sex"),
    "cand_det_sextractor_nnw": winter_file_dir.joinpath("default.nnw"),
    "cand_det_sextractor_filter": winter_file_dir.joinpath("default.conv"),
    "cand_det_sextractor_params": winter_file_dir.joinpath("Scorr.param"),
}

swarp_config_path = winter_file_dir.joinpath("config.swarp")
scamp_config_path = winter_file_dir.joinpath("astrom.scamp")
winter_mask_path = winter_file_dir.joinpath("winter_mask.fits")

winter_candidate_config = winter_file_dir.joinpath("candidates.sql")

psfex_path = winter_file_dir.joinpath("photom.psfex")

winter_cal_requirements = [
    CalRequirement(
        target_name="dark", required_field="EXPTIME", required_values=["120.0"]
    ),
]
