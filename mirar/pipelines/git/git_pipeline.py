"""
Module to run the SEDMv2 data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.git.blocks import (
    build_log,
    imsub,
    load_raw,
    load_stack_lt,
    reduce_raw_lt,
)
from mirar.pipelines.git.config import PIPELINE_NAME

logger = logging.getLogger(__name__)


class GITPipeline(Pipeline):
    """
    Class to run GIT/LT data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 30000  # no idea, for pylint
    all_pipeline_configurations = {
        "default": load_raw + build_log + imsub,
        "lt": reduce_raw_lt,
        "lt_sub": load_stack_lt + imsub,
    }

    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError
