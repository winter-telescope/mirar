"""
Module with typed mappings of the sextractor configuration files
"""

from pathlib import Path
from typing import TypedDict


class SextractorConfig(TypedDict):
    """Config file paths expanded into :class:`~.Sextractor`"""

    config_path: Path
    parameter_path: Path
    filter_path: Path
    starnnw_path: Path


class SextractorCandidateConfig(TypedDict):
    """Config file paths expanded into :class:`~.ZOGYSourceDetector`"""

    cand_det_sextractor_config: Path
    cand_det_sextractor_filter: Path
    cand_det_sextractor_nnw: Path
    cand_det_sextractor_params: Path
