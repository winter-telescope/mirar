"""
Module with typed mappings of the sextractor configuration files
"""

from pathlib import Path
from typing import TypedDict


class SextractorConfig(TypedDict):
    """Config file paths expanded into :class:`~.Sextractor`"""

    config_path: str | Path
    parameter_path: str | Path
    filter_path: str | Path
    starnnw_path: str | Path


class SextractorCandidateConfig(TypedDict):
    """Config file paths expanded into :class:`~.ZOGYSourceDetector`"""

    cand_det_sextractor_config: str | Path
    cand_det_sextractor_filter: str | Path
    cand_det_sextractor_nnw: str | Path
    cand_det_sextractor_params: str | Path
