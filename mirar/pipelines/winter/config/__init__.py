"""
This file contains the configuration for the winter pipeline.
"""

from pathlib import Path

from fastavro.schema import load_schema

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

sextractor_photometry_psf_config = {
    "config_path": winter_file_dir.joinpath("photomCat.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("photom_psf.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}


sextractor_anet_config = {
    "config_path_boardid_1_5": winter_file_dir.joinpath("astrom_anet_1_5.sex"),
    "config_path_boardid_0_2_3_4": winter_file_dir.joinpath("astrom_anet_0_2_3_4.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("astrom.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}

sextractor_reference_config = {
    "config_path": winter_file_dir.joinpath("reference.sex"),
    "parameter_path": winter_file_dir.joinpath("photom.param"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}

sextractor_candidate_config = {
    "cand_det_sextractor_config": winter_file_dir.joinpath("candidate_detection.sex"),
    "cand_det_sextractor_nnw": winter_file_dir.joinpath("default.nnw"),
    "cand_det_sextractor_filter": winter_file_dir.joinpath("default.conv"),
    "cand_det_sextractor_params": winter_file_dir.joinpath("Scorr.param"),
}

swarp_config_path = winter_file_dir.joinpath("config.swarp")
scamp_config_path = winter_file_dir.joinpath("astrom.scamp")
winter_mask_path = winter_file_dir.joinpath("winter_mask.fits")

psfex_path = winter_file_dir.joinpath("photom.psfex")

winter_cal_requirements = [
    CalRequirement(
        target_name="dark", required_field="EXPTIME", required_values=["120.0"]
    ),
]

winter_avro_schema_path = winter_file_dir.joinpath("avro_schema/winter.alert.avsc")
winter_avro_schema = load_schema(winter_avro_schema_path)
winter_prv_schema = winter_avro_schema["__named_schemas"]["winter.alert.prv_candidate"]
prv_candidate_cols = [x["name"] for x in winter_prv_schema["fields"]]
