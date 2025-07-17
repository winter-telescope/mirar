"""
Module to run the WASP data reduction pipeline
"""

import logging
from pathlib import Path

from mirar.data import Image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.lmi.blocks import build_log, load_raw, reduce, subtract
from mirar.pipelines.lmi.config import PIPELINE_NAME, lmi_cal_requirements
from mirar.pipelines.lmi.config.constants import LMI_NONLINEAR_LEVEL
from mirar.pipelines.lmi.load_lmi_image import load_raw_lmi_image

logger = logging.getLogger(__name__)


class LMIPipeline(Pipeline):
    """
    Class to run LMI data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = LMI_NONLINEAR_LEVEL  # no idea, for pylint
    all_pipeline_configurations = {
        "default": load_raw + reduce + subtract,
        "log": load_raw + build_log,
        "reduce": load_raw + reduce,
    }

    defalut_cal_requirements = lmi_cal_requirements

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError

    @staticmethod
    def _load_raw_image(path: str | Path) -> Image | list[Image]:
        return load_raw_lmi_image(path)
