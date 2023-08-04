"""
Module to run the SEDMv2 data reduction pipeline
"""
import logging
import os
from pathlib import Path

from mirar.data import Image
from mirar.downloader.caltech import download_via_ssh
from mirar.io import open_mef_image
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.sedmv2.blocks import (
    detrend_only,
    image_photometry,
    load_raw,
    process_stellar,
    process_transient,
)
from mirar.pipelines.sedmv2.config import PIPELINE_NAME, sedmv2_cal_requirements
from mirar.pipelines.sedmv2.load_sedmv2_image import load_raw_sedmv2_mef

sedmv2_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)


class SEDMv2Pipeline(Pipeline):
    """
    Class to run SEDMv2 data reduction pipeline
    """

    name = PIPELINE_NAME

    non_linear_level = 30000  # no idea, for pylint
    default_cal_requirements = sedmv2_cal_requirements
    all_pipeline_configurations = {
        "default": load_raw + process_stellar,
        "default_stellar": load_raw + process_stellar + image_photometry,
        "default_transient": load_raw + process_transient,  # +imsub,
        "realtime": load_raw + detrend_only,  # +much more...
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
    # def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
    #   return load_raw_sedmv2_image(path)
    def _load_raw_image(path: str | Path) -> Image | list[Image]:
        return open_mef_image(path, load_raw_sedmv2_mef)
