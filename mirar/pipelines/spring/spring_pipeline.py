"""
Module to run the SPRING data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.spring.blocks import (
load_raw
# TODO :  import the pipeline blocks HERE!!
)

from mirar.pipelines.spring.config import PIPELINE_NAME

logger = logging.getLogger(__name__)


class SPRINGPipeline(Pipeline):
    """
    Class to run GIT/LT data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 30000  # no idea, for pylint
    all_pipeline_configurations = {
        "default": [],
        "load_only": load_raw,
        "log": load_raw + csvlog,
    }

    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError
