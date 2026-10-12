"""
This file contains the configuration for the winter pipeline.
"""

import logging
import os
from pathlib import Path

from mirar.paths import base_raw_dir
from mirar.processors.astromatic.sextractor.config_type import (
    SextractorCandidateConfig,
    SextractorConfig,
)
from mirar.processors.skyportal.client import SkyportalClient
from mirar.processors.utils.cal_hunter import CalRequirement
from mirar.utils.avro_tools import load_avro_schema

logger = logging.getLogger(__name__)

PIPELINE_NAME = "winter"

winter_file_dir = Path(__file__).parent / "files"

sextractor_astrometry_config: SextractorConfig = {
    "config_path": winter_file_dir / "astrom.sex",
    "filter_path": winter_file_dir / "default.conv",
    "parameter_path": winter_file_dir / "astrom.param",
    "starnnw_path": winter_file_dir / "default.nnw",
}

sextractor_astromstats_config: SextractorConfig = {
    **sextractor_astrometry_config,
    "parameter_path": winter_file_dir / "astromstats.param",
}

sextractor_photometry_config: SextractorConfig = {
    "config_path": winter_file_dir / "photomCat.sex",
    "filter_path": winter_file_dir / "default.conv",
    "parameter_path": winter_file_dir / "photom.param",
    "starnnw_path": winter_file_dir / "default.nnw",
}

sextractor_photometry_psf_config: SextractorConfig = {
    "config_path": winter_file_dir / "photomCat.sex",
    "filter_path": winter_file_dir / "default.conv",
    "parameter_path": winter_file_dir / "photom_psf.param",
    "starnnw_path": winter_file_dir / "default.nnw",
}


sextractor_anet_config = {
    "config_path_boardid_1_5_6": winter_file_dir / "astrom_anet_1_5_6.sex",
    "config_path_boardid_2_3_4": winter_file_dir / "astrom_anet_2_3_4.sex",
    "filter_path": winter_file_dir / "default.conv",
    "parameter_path": winter_file_dir / "astrom.param",
    "starnnw_path": winter_file_dir / "default.nnw",
}

sextractor_reference_config: SextractorConfig = {
    "config_path": winter_file_dir / "reference.sex",
    "parameter_path": winter_file_dir / "photom.param",
    "filter_path": winter_file_dir / "default.conv",
    "starnnw_path": winter_file_dir / "default.nnw",
}

sextractor_reference_psf_phot_config: SextractorConfig = {
    "config_path": winter_file_dir / "reference.sex",
    "parameter_path": winter_file_dir / "photom_psf.param",
    "filter_path": winter_file_dir / "default.conv",
    "starnnw_path": winter_file_dir / "default.nnw",
}

sextractor_candidate_config: SextractorCandidateConfig = {
    "cand_det_sextractor_config": winter_file_dir / "candidate_detection.sex",
    "cand_det_sextractor_nnw": winter_file_dir / "default.nnw",
    "cand_det_sextractor_filter": winter_file_dir / "default.conv",
    "cand_det_sextractor_params": winter_file_dir / "Scorr.param",
}

swarp_config_path = winter_file_dir / "config.swarp"
scamp_config_path = winter_file_dir / "astrom.scamp"
winter_mask_path = winter_file_dir / "winter_mask.fits"

psfex_path = winter_file_dir / "photom.psfex"
ref_psfex_path = winter_file_dir / "reference.psfex"

winter_calhunter_directory = os.getenv("WINTER_CAL_DIRECTORY", None)
if winter_calhunter_directory is None:
    winter_calhunter_directory = base_raw_dir
    logger.debug(
        f"No specific calibration directory set. Using "
        f"{winter_calhunter_directory} as the directory for CalHunter. "
        f"If you want to change it, set the WINTER_CAL_DIRECTORY "
        f"environment variable."
    )
else:
    winter_calhunter_directory = Path(winter_calhunter_directory)

base_winter_cal_requirements = [
    CalRequirement(
        target_name="dark",
        required_field="EXPTIME",
        required_values=[
            "120.0",  # J/Y
            # "60.0",  # Hs
            # "3.0",  # J flats
            # "4.0",  # Y flats
            # "5.0",  # Hs flats
        ],
    ),
]

winter_cal_requirements = base_winter_cal_requirements

winter_avro_schema_path = winter_file_dir / "avro_schema/winter.alert.avsc"
winter_avro_schema = load_avro_schema(winter_avro_schema_path)
winter_prv_schema = winter_avro_schema["__named_schemas"]["winter.alert.prv_candidate"]
prv_candidate_cols = [x["name"] for x in winter_prv_schema["fields"]]

# To send to preview.fritz.science
winter_preview_config = {
    "origin": "mirar",
    "group_ids": [1092],
    "fritz_filter_id": 1018,
    "instrument_id": 1066,
    "stream_id": 1008,
    "update_thumbnails": True,
    "skyportal_client": SkyportalClient(base_url="https://preview.fritz.science/api/"),
}

# To send to fritz.science
winter_fritz_config = {
    "origin": "mirar",
    "group_ids": [1657],
    "fritz_filter_id": 1185,
    "instrument_id": 1087,
    "stream_id": 1005,
    "update_thumbnails": True,
    "skyportal_client": SkyportalClient(base_url="https://fritz.science/api/"),
    "annotation_keys": [
        "rb",
        "chipsf",
        "fwhm",
        "scorr",
        "nneg",
        "mindtoedge",
        "diffmaglim",
        "distpsnr1",
        "sgmag1",
        "srmag1",
        "simag1",
    ],
}
