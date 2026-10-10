"""
Module containing SPRING-specific paths
"""

import os
from pathlib import Path

from mirar.pipelines.spring.config.constants import (
    DB_NAME,
    PIPELINE_NAME,
    SPRING_GAIN,
    SPRING_PIXEL_SCALE,
)
from mirar.processors.astromatic.sextractor.config_type import (
    SextractorCandidateConfig,
    SextractorConfig,
)
from mirar.processors.skyportal.client import SkyportalClient
from mirar.processors.utils.cal_hunter import CalRequirement

spring_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(spring_dir, "files")

sextractor_photometry_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "photomCat.sex"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "parameter_path": Path(astromatic_config_dir, "photom.param"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}

sextractor_PSF_photometry_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "photomCat.sex"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "parameter_path": Path(astromatic_config_dir, "photomPSF.param"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}

sextractor_astrometry_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "astrom.sex"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "parameter_path": Path(astromatic_config_dir, "astrom.param"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}


sextractor_candidates_config: SextractorCandidateConfig = {
    "cand_det_sextractor_config": Path(astromatic_config_dir, "photomCat.sex"),
    "cand_det_sextractor_nnw": Path(astromatic_config_dir, "default.nnw"),
    "cand_det_sextractor_filter": Path(astromatic_config_dir, "default.conv"),
    "cand_det_sextractor_params": Path(astromatic_config_dir, "Scorr.param"),
}

sextractor_reference_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "reference.sex"),
    "parameter_path": Path(astromatic_config_dir, "photom.param"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}

sextractor_reference_psf_phot_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "reference.sex"),
    "parameter_path": Path(astromatic_config_dir, "photomPSF_ref.param"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}

scamp_path = os.path.join(astromatic_config_dir, "scamp.conf")

swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")

psfex_config_path = os.path.join(astromatic_config_dir, "photom.psfex")

psfex_sci_config_path = os.path.join(astromatic_config_dir, "photom_sci.psfex")

ref_psfex_path = os.path.join(astromatic_config_dir, "reference.psfex")

spring_cal_requirements = [
    CalRequirement(
        target_name="dark",
        required_field="EXPTIME",
        required_values=[
            "60.0",  # J
        ],
    )
]

spring_fritz_config = {
    "group_ids": [1],
    "instrument_id": 1138,
    "update_thumbnails": False,
    "skyportal_client": SkyportalClient(base_url="https://fritz.science/api/"),
}
