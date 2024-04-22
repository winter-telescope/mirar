"""
This file contains the configuration for the winter pipeline.
"""

from pathlib import Path

from fastavro.schema import load_schema

from mirar.processors.utils.cal_hunter import CalRequirement

PIPELINE_NAME = "nires"

file_dir = Path(__file__).parent.joinpath("files")

sextractor_astrometry_config = {
    "config_path": file_dir.joinpath("astrom.sex"),
    "filter_path": file_dir.joinpath("default.conv"),
    "parameter_path": file_dir.joinpath("astrom.param"),
    "starnnw_path": file_dir.joinpath("default.nnw"),
}
