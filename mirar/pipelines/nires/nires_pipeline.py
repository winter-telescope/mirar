"""
Module to run the NIRES acquisition camera data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.nires.blocks import build_log, flat_calibrate, load_raw

PIPELINE_NAME = "nires"

logger = logging.getLogger(__name__)


class NIRESPipeline(Pipeline):
    """
    Class to run NIRES data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 30000  # no idea, for pylint
    all_pipeline_configurations = {"default": load_raw + build_log + flat_calibrate}

    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError
