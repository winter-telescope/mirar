"""
Module with pipeline for building reference images in the IR from WFAU
"""
import logging

from mirar.downloader.caltech import download_via_ssh
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.winter.blocks import (
    commissioning,
    commissioning_dark,
    commissioning_flat,
    commissioning_multiboard_stack,
    commissioning_noise,
    commissioning_proc,
    commissioning_reduce,
    commissioning_stack,
    full_commissioning,
    log,
    refbuild,
)
from mirar.pipelines.winter.config import PIPELINE_NAME

logger = logging.getLogger(__name__)


class WINTERPipeline(Pipeline):
    """
    Pipeline for processing WINTER data
    """

    name = "winter"

    all_pipeline_configurations = {
        "refbuild": refbuild,
        "commissioning": commissioning,
        "log": log,
        "commissioning_dark": commissioning_dark,
        "commissioning_stack": commissioning_stack,
        "commissioning_flat": commissioning_flat,
        "commissioning_reduce": commissioning_reduce,
        "commissioning_proc": commissioning_proc,
        "commissioning_multiboard_stack": commissioning_multiboard_stack,
        "full_commissioning": full_commissioning,
        "commissioning_noise": commissioning_noise,
    }

    gain = 1.0
    non_linear_level = 65535

    @staticmethod
    def _load_raw_image(path: str):
        pass

    @staticmethod
    def download_raw_images_for_night(night: str):
        download_via_ssh(
            server="winter.caltech.edu",
            base_dir="/data/loki/raw_data/winter",
            night=night,
            pipeline=PIPELINE_NAME,
            server_sub_dir="raw",
        )
