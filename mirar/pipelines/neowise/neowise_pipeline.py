"""
Module to run the NEOWISE data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.neowise.blocks import (
    image_photometry,
    load_ref_neowise,
    load_stack_neowise,
    save_image,
    subtract,
)
from mirar.pipelines.neowise.config import PIPELINE_NAME

logger = logging.getLogger(__name__)


class NEOWISEPipeline(Pipeline):
    """
    Class to run NEOWISE image subtraction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 30000  # no idea, for pylint
    all_pipeline_configurations = {
        "subtract": load_stack_neowise + subtract + image_photometry,
        "prepare": load_ref_neowise + save_image,
    }

    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError
