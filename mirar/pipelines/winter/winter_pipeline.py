"""
Module with pipeline for building reference images in the IR from WFAU
"""
import logging

from mirar.data import Image
from mirar.downloader.caltech import download_via_ssh
from mirar.io import open_mef_image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.winter.blocks import (
    build_test,
    csvlog,
    detect_candidates,
    detrend_unpacked,
    imsub,
    load_candidates,
    load_stack,
    load_test,
    load_test_stack,
    only_ref,
    photcal_stacks,
    process_candidates,
    realtime,
    reduce,
    reduce_unpacked,
    refbuild,
    reftest,
    unpack_all,
    unpack_subset,
)
from mirar.pipelines.winter.config import PIPELINE_NAME, winter_cal_requirements
from mirar.pipelines.winter.load_winter_image import load_raw_winter_mef

logger = logging.getLogger(__name__)


class WINTERPipeline(Pipeline):
    """
    Pipeline for processing WINTER data
    """

    name = "winter"
    default_cal_requirements = winter_cal_requirements

    all_pipeline_configurations = {
        "unpack_subset": unpack_subset,
        "unpack_all": unpack_all,
        "detrend_unpacked": detrend_unpacked,
        "reduce": reduce,
        "reduce_unpacked": reduce_unpacked,
        "photcal_stacks": photcal_stacks,
        "buildtest": build_test,
        "test": load_test + csvlog + realtime,
        "test_imsub": load_test_stack + imsub + detect_candidates + process_candidates,
        "refbuild": refbuild,
        "reftest": reftest,
        "only_ref": only_ref,
        "realtime": realtime,
        "imsub": load_stack + imsub,
        "detect_candidates": load_stack + imsub + detect_candidates,
        "process_candidates": load_candidates + process_candidates,
        "full_imsub": load_stack + imsub + detect_candidates + process_candidates,
    }

    non_linear_level = 40000.0

    @staticmethod
    def _load_raw_image(path: str) -> Image | list[Image]:
        return open_mef_image(path, load_raw_winter_mef, extension_key="BOARD_ID")

    @staticmethod
    def download_raw_images_for_night(night: str):
        download_via_ssh(
            server="winter.caltech.edu",
            base_dir="/data/loki/raw_data/winter",
            night=night,
            pipeline=PIPELINE_NAME,
            server_sub_dir="raw",
        )
