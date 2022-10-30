import astropy.io.fits
import numpy as np
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.downloader.caltech import download_via_ssh
from winterdrp.pipelines.wirc.blocks import load_raw, reduce, imsub
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
import logging

logger = logging.getLogger(__name__)


pipeline_name = "wirc"


class WircPipeline(Pipeline):

    name = pipeline_name

    non_linear_level = 30000
    gain = 1.2

    # Set up elements to use

    header_keys = [
        "UTSHUT",
        'OBJECT',
        "FILTER",
        "EXPTIME",
        "COADDS",
    ]

    all_pipeline_configurations = {
        "default": load_raw + reduce,
        "imsub": load_raw + imsub
    }

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        download_via_ssh(
            server="gayatri.caltech.edu",
            base_dir="/scr2/ptf/observation_data",
            prefix="WIRC_",
            night=night,
            pipeline=pipeline_name,
            server_sub_dir="raw"
        )

    @staticmethod
    def load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_wirc_image(path)

