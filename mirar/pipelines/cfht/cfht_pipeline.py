"""
Module to run the SEDMv2 data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.cfht.blocks import astrometry
from mirar.pipelines.cfht.config import PIPELINE_NAME

logger = logging.getLogger(__name__)


class CFHTPipeline(Pipeline):
    """
    Class to run GIT/LT data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 30000  # no idea, for pylint
    all_pipeline_configurations = {
        "default": astrometry,
    }

    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError
