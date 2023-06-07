"""
Module to run the summer data reduction pipeline
"""
import logging
import os

import astropy.io.fits
import numpy as np

from mirar.downloader.caltech import download_via_ssh
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.summer.blocks import (
    build_log,
    cal_hunter,
    export_raw,
    imsub,
    load_processed,
    load_raw,
    load_test,
    load_test_proc,
    process_raw,
    sim_realtime,
    subtract,
    test_cr,
)
from mirar.pipelines.summer.config import PIPELINE_NAME, summer_cal_requirements
from mirar.pipelines.summer.load_summer_image import load_raw_summer_image

summer_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)


class SummerPipeline(Pipeline):
    """
    Class to run summer data reduction pipeline
    """

    name = PIPELINE_NAME
    default_cal_requirements = summer_cal_requirements

    all_pipeline_configurations = {
        "default": load_raw + build_log + export_raw + cal_hunter + process_raw,
        "test": load_test + export_raw + process_raw,
        "postprocess": build_log,
        "imsub": load_processed + imsub,
        "test_imsub": load_test_proc + subtract,
        "full": load_raw + build_log + export_raw + cal_hunter + process_raw + imsub,
        "realtime": export_raw + process_raw,
        "log": load_raw + build_log,
        "simrealtime": sim_realtime,
        "testlog": load_test + build_log,
        "crtest": load_raw + test_cr + build_log,
        "dbtest": load_raw + export_raw,
    }

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        download_via_ssh(
            server="winter.caltech.edu",
            base_dir="/data/loki/raw_data/summer",
            night=night,
            pipeline=PIPELINE_NAME,
            server_sub_dir="raw",
        )

    @staticmethod
    def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_summer_image(path)
