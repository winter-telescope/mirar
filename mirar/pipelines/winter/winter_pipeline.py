"""
Module with pipeline for building reference images in the IR from WFAU
"""
import logging

from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.winter.blocks import refbuild

logger = logging.getLogger(__name__)


class WINTERPipeline(Pipeline):
    """
    Pipeline for building reference images in the IR from WFAU
    """

    name = "winter"

    all_pipeline_configurations = {"refbuild": refbuild}

    gain = 1.0
    non_linear_level = 65535

    @staticmethod
    def _load_raw_image(path: str):
        pass

    @staticmethod
    def download_raw_images_for_night(night: str):
        pass
