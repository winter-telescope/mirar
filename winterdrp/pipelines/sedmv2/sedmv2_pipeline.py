"""
Module to run the SEDMv2 data reduction pipeline
"""
import logging
import os

import astropy.io.fits
import numpy as np

from winterdrp.downloader.caltech import download_via_ssh
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.pipelines.sedmv2.blocks import (
    build_log,
    cal_hunter,
    imsub,
    load_raw,
    load_test,
    process_raw,
    sim_realtime,
)
from winterdrp.pipelines.sedmv2.config import PIPELINE_NAME, sedmv2_cal_requirements
from winterdrp.pipelines.sedmv2.load_sedmv2_image import load_raw_sedmv2_image

sedmv2_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)


class SEDMv2Pipeline(Pipeline):
    """
    Class to run SEDMv2 data reduction pipeline
    """

    name = PIPELINE_NAME
    gain = 1
    non_linear_level = 30000  # no idea, for pylint
    default_cal_requirements = sedmv2_cal_requirements
    # removed export_raw and load_processed blocks
    all_pipeline_configurations = {
        "default": load_raw + cal_hunter + process_raw,
        "test": load_test + process_raw,
        "postprocess": build_log,
        "imsub": imsub,
        "full": load_raw + build_log + cal_hunter + process_raw + imsub,
        "realtime": process_raw,
        "log": load_raw + build_log,
        "simrealtime": sim_realtime,
        "testlog": load_test + build_log,
    }

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        download_via_ssh(
            server="jagati.caltech.edu",
            base_dir="/data/viraj/winter_data/commissioning/raw/",
            night=night,
            pipeline=PIPELINE_NAME,
        )

    @staticmethod
    def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_sedmv2_image(path)
