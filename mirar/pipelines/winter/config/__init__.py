from pathlib import Path

PIPELINE_NAME = "winter"

winter_file_dir = Path(__file__).parent.joinpath("files")

sextractor_astrometry_config = {
    "config_path": winter_file_dir.joinpath("photomCat.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("astrom.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}


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

sextractor_autoastrometry_config = {
    "config_path": winter_file_dir.joinpath("astrom.sex"),
    "filter_path": winter_file_dir.joinpath("default.conv"),
    "parameter_path": winter_file_dir.joinpath("astrom.param"),
    "starnnw_path": winter_file_dir.joinpath("default.nnw"),
}

swarp_config_path = winter_file_dir.joinpath("swarp.config")
