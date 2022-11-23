import logging
import os
import astropy.io.fits
import numpy as np

from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.downloader.caltech import download_via_ssh

from winterdrp.pipelines.summer.config import PIPELINE_NAME, summer_cal_requirements
from winterdrp.pipelines.summer.load_summer_image import load_raw_summer_image
from winterdrp.pipelines.summer.blocks import load_raw, build_log, load_processed, imsub, export_raw, \
    cal_hunter, process_raw, load_test, load_test_proc, subtract, sim_realtime

summer_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)


class SummerPipeline(Pipeline):

    name = PIPELINE_NAME
    default_cal_requirements = summer_cal_requirements

    all_pipeline_configurations = {
        "default": load_raw + export_raw + cal_hunter + process_raw,
        "test": load_test + export_raw + process_raw,
        "postprocess": build_log,
        'imsub': load_processed + imsub,
        "test_imsub": load_test_proc + subtract,
        "full": load_raw + build_log + export_raw + cal_hunter + process_raw + imsub,
        "realtime": export_raw + process_raw,
        "log": load_raw + build_log,
        "simrealtime": sim_realtime,
        "testlog": load_test + build_log
    }

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        download_via_ssh(
            server="jagati.caltech.edu",
            base_dir="/data/viraj/winter_data/commissioning/raw/",
            night=night,
            pipeline=PIPELINE_NAME
        )

    @staticmethod
    def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_summer_image(path)
