"""
Module with pipeline for building reference images in the IR from WFAU
"""

import logging

from mirar.data import Image
from mirar.downloader.caltech import download_via_ssh
from mirar.io import open_mef_image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.winter.blocks import (
    astrometry,
    avro_broadcast,
    avro_export,
    build_test,
    csvlog,
    detect_candidates,
    detrend_unpacked,
    diff_forced_photometry,
    extract_all,
    focus_cals,
    full_reduction,
    imsub,
    load_avro,
    load_calibrated,
    load_final_stack,
    load_raw,
    load_skyportal,
    load_test,
    mask_and_split,
    mosaic,
    only_ref,
    photcal_stacks,
    plot_stack,
    process_candidates,
    realtime,
    reduce,
    reduce_unpacked,
    reduce_unpacked_subset,
    refbuild,
    reftest,
    save_raw,
    select_split_subset,
    send_to_skyportal,
    stack_forced_photometry,
    stack_stacks,
    unpack_all,
    unpack_subset,
)
from mirar.pipelines.winter.config import PIPELINE_NAME, winter_cal_requirements
from mirar.pipelines.winter.load_winter_image import load_raw_winter_mef
from mirar.pipelines.winter.models import set_up_winter_databases

logger = logging.getLogger(__name__)


class WINTERPipeline(Pipeline):
    """
    Pipeline for processing WINTER data
    """

    name = "winter"
    default_cal_requirements = winter_cal_requirements

    all_pipeline_configurations = {
        "astrometry": load_calibrated + astrometry,
        "unpack_subset": unpack_subset,
        "unpack_all": unpack_all,
        "detrend_unpacked": detrend_unpacked,
        "imsub": load_final_stack + imsub,
        "reduce": reduce,
        "reduce_unpacked": reduce_unpacked,
        "photcal_stacks": photcal_stacks,
        "plot_stacks": load_final_stack + plot_stack,
        "buildtest": build_test,
        "test": load_test
        + csvlog
        + extract_all
        + mask_and_split
        + select_split_subset
        + save_raw
        + full_reduction
        + imsub
        + detect_candidates
        + process_candidates,
        "refbuild": refbuild,
        "reftest": reftest,
        "only_ref": only_ref,
        "realtime": realtime,
        "detect_candidates": load_final_stack
        + imsub
        + detect_candidates
        + process_candidates
        + avro_broadcast,
        "default": reduce
        + imsub
        + detect_candidates
        + process_candidates
        + avro_broadcast,
        "default_subset": reduce_unpacked_subset
        + imsub
        + detect_candidates
        + process_candidates
        + avro_broadcast,
        "stack_stacks": load_final_stack + stack_stacks,
        "stack_stacks_db": stack_stacks,
        "focus_cals": focus_cals,
        "mosaic": mosaic,
        "log": load_raw + extract_all + csvlog,
        "skyportal": load_skyportal + send_to_skyportal,
        "diff_forced_phot": diff_forced_photometry,
        "stack_forced_phot": stack_forced_photometry,
        "rebroadcast_avro": load_avro + avro_export,
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

    def set_up_pipeline(self):
        set_up_winter_databases()
