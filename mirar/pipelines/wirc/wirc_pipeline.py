"""
Module for the WIRC (https://doi.org/10.1117/12.460336) pipeline
"""

import logging
from pathlib import Path

from mirar.data import Image
from mirar.downloader.caltech import download_via_ssh
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.wirc.blocks import (
    candidates,
    imsub,
    load_raw,
    load_stack,
    log,
    reduce,
    write_file,
    load_mosfire_raw,
    load_hst_raw,
    psfex,
)
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image

logger = logging.getLogger(__name__)


PIPELINE_NAME = "wirc"


class WircPipeline(Pipeline):
    """
    Pipeline for WIRC on the Palomar 200-inch telescope
    """

    name = PIPELINE_NAME

    non_linear_level = 30000

    all_pipeline_configurations = {
        "default": load_raw + reduce + imsub,
        "reduce": load_raw + reduce,
        "imsub": load_stack + imsub,
        "imsub_mosfire": load_mosfire_raw + imsub,
        "log": load_raw + log,
        "candidates": candidates,
        "pass_through_pipe": load_mosfire_raw + write_file,
        "pass_hst_through_pipe": load_hst_raw + write_file,
        "psfex":load_mosfire_raw + psfex,
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
    def _load_raw_image(path: str | Path) -> Image | list[Image]:
        return load_raw_wirc_image(path)
