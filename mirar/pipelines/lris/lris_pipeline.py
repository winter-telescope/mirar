"""
Module to run the WASP data reduction pipeline
"""

import logging
from pathlib import Path

from mirar.data import Image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.lris.blocks import build_log, load_raw, reduce, subtract
from mirar.pipelines.lris.config import PIPELINE_NAME, lris_cal_requirements
from mirar.pipelines.lris.config.constants import LRIS_NONLINEAR_LEVEL
from mirar.pipelines.lris.load_lris_image import load_raw_lris_image

logger = logging.getLogger(__name__)


class LRISPipeline(Pipeline):
    """
    Class to run LRIS data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = LRIS_NONLINEAR_LEVEL  # no idea, for pylint
    all_pipeline_configurations = {
        "default": load_raw + reduce + subtract,
        "log": load_raw + build_log,
        "reduce": load_raw + reduce,
    }

    defalut_cal_requirements = lris_cal_requirements

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError

    @staticmethod
    def _load_raw_image(path: str | Path) -> Image | list[Image]:
        return load_raw_lris_image(path)
