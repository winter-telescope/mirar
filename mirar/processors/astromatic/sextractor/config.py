"""
Config Type for Sextractor Processor
"""

from pathlib import Path
from typing import TypedDict


class SextractorConfig(TypedDict):
    config_path: str | Path
    filter_path: str | Path
    parameter_path: str | Path
    starnnw_path: str | Path
