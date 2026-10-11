"""
Module with typed mappings of the sextractor configuration files
"""

from pathlib import Path
from typing import TypedDict


class SextractorConfig(TypedDict):
    """
    Config file paths expanded into :class:`~.Sextractor`

    :param config_path: Path of the sextractor config file
    :param parameter_path: Path of the sextractor parameter file
    :param filter_path: Path of the sextractor filter file
    :param starnnw_path: Path of the sextractor starnnw file
    """

    config_path: Path
    parameter_path: Path
    filter_path: Path
    starnnw_path: Path


class SextractorCandidateConfig(TypedDict):
    """
    Config file paths expanded into :class:`~.ZOGYSourceDetector`

    :param cand_det_sextractor_config: Path of the sextractor config file
    :param cand_det_sextractor_filter: Path of the sextractor filter file
    :param cand_det_sextractor_nnw: Path of the sextractor starnnw file
    :param cand_det_sextractor_params: Path of the sextractor parameter file
    """

    cand_det_sextractor_config: Path
    cand_det_sextractor_filter: Path
    cand_det_sextractor_nnw: Path
    cand_det_sextractor_params: Path
