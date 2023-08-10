"""
Module containing WIRC-specific paths
"""
from pathlib import Path

from fastavro.schema import load_schema

wirc_file_dir = Path(__file__).parent.joinpath("files")
wirc_schema_dir = Path(__file__).parent.joinpath("schema")

wirc_mask_path = wirc_file_dir.joinpath("wirc_bad_feb_2018.fits")

sextractor_astrometry_config = {
    "config_path": wirc_file_dir.joinpath("stack.sex"),
    "filter_path": wirc_file_dir.joinpath("default.conv"),
    "parameter_path": wirc_file_dir.joinpath("astrom.param"),
    "starnnw_path": wirc_file_dir.joinpath("default.nnw"),
}

sextractor_photometry_config = {
    "config_path": wirc_file_dir.joinpath("matchcat.sex"),
    "filter_path": wirc_file_dir.joinpath("default.conv"),
    "parameter_path": wirc_file_dir.joinpath("astrom.param"),
    "starnnw_path": wirc_file_dir.joinpath("default.nnw"),
}

scamp_fp_path = wirc_file_dir.joinpath("scamp_fp.conf")

swarp_sp_path = wirc_file_dir.joinpath("second_pass.swarp")

sextractor_reference_config = {
    "config_path": wirc_file_dir.joinpath("photomCat.sex"),
    "parameter_path": wirc_file_dir.joinpath("photom.param"),
    "filter_path": wirc_file_dir.joinpath("default.conv"),
    "starnnw_path": wirc_file_dir.joinpath("default.nnw"),
}

sextractor_candidate_config = {
    "cand_det_sextractor_config": wirc_file_dir.joinpath("photomCat.sex"),
    "cand_det_sextractor_nnw": wirc_file_dir.joinpath("default.nnw"),
    "cand_det_sextractor_filter": wirc_file_dir.joinpath("default.conv"),
    "cand_det_sextractor_params": wirc_file_dir.joinpath("Scorr.param"),
}

psfex_path = wirc_file_dir.joinpath("photom.psfex")

wirc_avro_schema_path = Path(__file__).parent.joinpath("avro_schema/wirc.alert.avsc")
wirc_avro_schema = load_schema(wirc_avro_schema_path)
wirc_prv_schema = wirc_avro_schema["__named_schemas"]["wirc.alert.prv_candidate"]
prv_candidate_cols = [x["name"] for x in wirc_prv_schema["fields"]]
