"""
Module with pipeline for building reference images in the IR from WFAU
"""
import logging

from mirar.data import Image
from mirar.downloader.caltech import download_via_ssh
from mirar.io import open_mef_image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.winter.blocks import (
    commissioning_multiboard_stack,
    commissioning_photcal,
    commissioning_photcal_indiv,
    final,
    full_commissioning,
    full_commissioning_all_boards,
    imsub,
    only_ref,
    realtime,
    reduce,
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
        "refbuild": refbuild,
        "commissioning_multiboard_stack": commissioning_multiboard_stack,
        "full_commissioning": full_commissioning,
        "commissioning_photcal": commissioning_photcal,
        "commissioning_photcal_indiv": commissioning_photcal_indiv,
        "full_commissioning_all_boards": full_commissioning_all_boards,
        "unpack_subset": unpack_subset,
        "unpack_all": unpack_all,
        "imsub": imsub,
        "final": final,
        "reduce": reduce,
        "reftest": reftest,
        "only_ref": only_ref,
        "realtime": realtime,
    }

    gain = 1.0
    non_linear_level = 65535

    @staticmethod
    def _load_raw_image(path: str) -> Image | list[Image]:
        return open_mef_image(path, load_raw_winter_mef)

    @staticmethod
    def download_raw_images_for_night(night: str):
        download_via_ssh(
            server="winter.caltech.edu",
            base_dir="/data/loki/raw_data/winter",
            night=night,
            pipeline=PIPELINE_NAME,
            server_sub_dir="raw",
        )
