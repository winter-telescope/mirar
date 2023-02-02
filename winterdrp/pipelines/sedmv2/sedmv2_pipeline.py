import logging
import os
import astropy.io.fits
import numpy as np

from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.downloader.caltech import download_via_ssh

from winterdrp.pipelines.sedmv2.config import PIPELINE_NAME, sedmv2_cal_requirements
from winterdrp.pipelines.sedmv2.load_sedmv2_image import load_raw_sedmv2_image
from winterdrp.pipelines.sedmv2.blocks import load_raw, build_log, imsub,\
    cal_hunter, process_raw, load_test, load_test_proc, subtract, sim_realtime

sedmv2_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)


class SEDMv2Pipeline(Pipeline):

    name = PIPELINE_NAME
    default_cal_requirements = sedmv2_cal_requirements
    # removed export_raw and load_processed blocks
    all_pipeline_configurations = {
        "default": load_raw + cal_hunter + process_raw,
        "test": load_test + process_raw,
        "postprocess": build_log,
        'imsub': imsub,
        "test_imsub": load_test_proc + subtract,
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
            base_dir="/data/viraj/winter_data/commissioning/raw/", # TODO: change directory/server
            night=night,
            pipeline=PIPELINE_NAME
        )

    @staticmethod
    def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_sedmv2_image(path)
