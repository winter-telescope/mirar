"""
Module to run the SPRING data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.spring.blocks import (
    astrometry,
    csvlog,
    dark_calibrate,
    flat_calibrate,
    load_raw,
    photcal_with_color,
    photcal_without_color,
    stack_dithers,
)
from mirar.pipelines.spring.config.constants import PIPELINE_NAME
from mirar.pipelines.spring.models import set_up_spring_databases

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
        "darkcal": load_raw + csvlog + dark_calibrate,
        "darks_flats": load_raw + csvlog + dark_calibrate + flat_calibrate,
        "astrometry": load_raw + dark_calibrate + flat_calibrate + astrometry,
        "stacking": load_raw
        + dark_calibrate
        + flat_calibrate
        + astrometry
        + stack_dithers,
        "photometry": load_raw
        + dark_calibrate
        + flat_calibrate
        + astrometry
        + stack_dithers
        + photcal_without_color,
        "photometry_color": load_raw
        + dark_calibrate
        + flat_calibrate
        + astrometry
        + stack_dithers
        + photcal_with_color,
    }

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError

    def set_up_pipeline(self):
        set_up_spring_databases()
