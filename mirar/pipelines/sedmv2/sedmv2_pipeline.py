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
from mirar.pipelines.sedmv2.blocks import (  # transient_phot_psfexsex,
    image_photometry,
    imsub,
    load_raw,
    process_all_psf_then_cal,
    process_all_psf_then_cal_no_color_term,
    process_stellar,
    process_transient_using_magauto,
    psf_phot_after_imsub,
    reduce_not0,
    transient_phot,
    upload_fritz,
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
        "default_transient": load_raw
        + process_transient_using_magauto
        + transient_phot,  # +imsub,
        "realtime": load_raw + reduce_not0,  # +much more...
        "psf_all_then_calibrate": load_raw + process_all_psf_then_cal,
        "transient_upload": load_raw
        + process_transient_using_magauto
        + transient_phot
        + upload_fritz,
        "transient_PSF": load_raw
        + process_all_psf_then_cal
        + imsub
        + psf_phot_after_imsub,  # + upload_fritz
        "transient_PSF_nocolor": load_raw
        + process_all_psf_then_cal_no_color_term
        + imsub
        + psf_phot_after_imsub,
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
