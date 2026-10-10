"""
Module containing WASP-specific paths
"""

# pylint: disable=duplicate-code

from pathlib import Path

from mirar.pipelines.wasp.config.constants import PIPELINE_NAME
from mirar.processors.astromatic.sextractor.config_type import SextractorConfig
from mirar.processors.utils.cal_hunter import CalRequirement

wasp_dir = Path(__file__).parent

astromatic_config_dir = wasp_dir / "files"

sextractor_photometry_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "photomCat.sex"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "parameter_path": Path(astromatic_config_dir, "photom.param"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}

sextractor_astrometry_config: SextractorConfig = {
    "config_path": Path(astromatic_config_dir, "astrom.sex"),
    "filter_path": Path(astromatic_config_dir, "default.conv"),
    "parameter_path": Path(astromatic_config_dir, "astrom.param"),
    "starnnw_path": Path(astromatic_config_dir, "default.nnw"),
}

scamp_path = astromatic_config_dir / "scamp.conf"

swarp_config_path = astromatic_config_dir / "config.swarp"

psfex_sci_config_path = astromatic_config_dir / "photom_sci.psfex"

wasp_cal_requirements = [
    CalRequirement(
        target_name="bias", required_field="EXPTIME", required_values=["0.0"]
    ),
]
