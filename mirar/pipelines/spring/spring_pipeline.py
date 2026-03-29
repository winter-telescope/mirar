"""
Module to run the SPRING data reduction pipeline
"""

import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.spring.blocks import (
    astrometry,
    csvlog,
    dark_calibrate,
    detect_candidates,
    diff_forced_photometry,
    flat_calibrate,
    imsub,
    load_final_stack,
    load_raw,
    load_sub,
    photcal_and_export,
    photcal_color_and_export,
    photcal_with_color,
    photcal_without_color,
    reduce,
    stack_dithers,
    stack_forced_photometry,
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
        "flats": load_raw + csvlog + flat_calibrate,
        "darks_flats": reduce,
        "astrometry": reduce + astrometry,
        "stacking": reduce + astrometry + stack_dithers,
        "photometry": reduce + astrometry + stack_dithers + photcal_without_color,
        "photometry_color": reduce + astrometry + stack_dithers + photcal_with_color,
        "photcal_and_export": photcal_and_export,
        "photcal_color_and_export": photcal_color_and_export,
        "subtraction": load_final_stack + imsub,
        "full_stack_fp": reduce + astrometry + stack_dithers + stack_forced_photometry,
        "full_diff_fp": reduce
        + astrometry
        + stack_dithers
        + imsub
        + diff_forced_photometry,
        "fp_from_stack": load_final_stack + stack_forced_photometry,
        "fp_from_diff": load_sub + diff_forced_photometry,
        "candidates_from_diff": load_sub + detect_candidates,
    }

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        """
        Download raw images for a night
        """
        raise NotImplementedError

    def set_up_pipeline(self):
        set_up_spring_databases()
