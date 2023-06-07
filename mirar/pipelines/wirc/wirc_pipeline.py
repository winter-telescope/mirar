"""
Module for the WIRC (https://doi.org/10.1117/12.460336) pipeline
"""
import logging

import astropy.io.fits
import numpy as np

from mirar.downloader.caltech import download_via_ssh
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.wirc.blocks import imsub, load_raw, reduce
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image

logger = logging.getLogger(__name__)


PIPELINE_NAME = "wirc"


class WircPipeline(Pipeline):
    """
    Pipeline for WIRC on the Palomar 200-inch telescope
    """

    name = PIPELINE_NAME

    non_linear_level = 30000
    gain = 1.2

    all_pipeline_configurations = {
        "default": load_raw + reduce,
        "imsub": load_raw + imsub,
    }

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        download_via_ssh(
            server="gayatri.caltech.edu",
            base_dir="/scr2/ptf/observation_data",
            prefix="WIRC_",
            night=night,
            pipeline=PIPELINE_NAME,
            server_sub_dir="raw",
        )

    @staticmethod
    def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_wirc_image(path)
