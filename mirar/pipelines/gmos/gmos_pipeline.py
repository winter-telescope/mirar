"""
Module to run the GMOS data reduction pipeline
"""

import logging
from pathlib import Path

from mirar.data import Image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.gmos.blocks import build_log, load_raw, reduce, subtract
from mirar.pipelines.gmos.config import PIPELINE_NAME
from mirar.pipelines.gmos.load_gmos_image import load_detrended_gmos_image

logger = logging.getLogger(__name__)


class GMOSPipeline(Pipeline):
    """
    Class to run GMOS data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 100000  # no idea, for pylint
    all_pipeline_configurations = {
        "default": load_raw + reduce + subtract,
        "log": load_raw + build_log,
        "reduce": load_raw + reduce,
    }

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError

    @staticmethod
    def _load_raw_image(path: str | Path) -> Image | list[Image]:
        return load_detrended_gmos_image(path)
