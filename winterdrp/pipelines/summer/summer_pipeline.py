import os
import astropy.io.fits
import numpy as np
from astropy.io.fits import HDUList
from winterdrp.pipelines.base_pipeline import Pipeline

from winterdrp.processors.bias import BiasCalibrator
from winterdrp.processors.flat import FlatCalibrator
from winterdrp.processors.utils import ImageSaver
from winterdrp.processors.astromatic import SextractorRunner

from winterdrp.pipelines.summer.calibration import select_bias, select_flats_archival

summer_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class SummerPipeline(Pipeline):

    name = "summer"

    astrometry_cal = ("GAIA", 9., 13.)
    photometry_cal = {
        "J": ()
    }

    # Set up elements to use

    header_keys = [
        "UTC",
        'FIELDID',
        "FILTERID",
        "EXPTIME",
        "OBSTYPE"
    ]

    batch_split_keys = ["RAWIMAGEPATH"]

    pipeline_configurations = {
        None: [
            (BiasCalibrator, select_bias),
            (FlatCalibrator, select_flats_archival),
            (ImageSaver, "preprocess"),
            (SextractorRunner, "pass1"),
            # "stack",
            # "dither"
        ]
    }

    @staticmethod
    def reformat_raw_data(
            img: HDUList,
            path: str
    ) -> [np.array, astropy.io.fits.Header]:
        header = img[0].header
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "SCIENCE"]
        header["BASENAME"] = os.path.basename(path)
        header["CALSTEPS"] = ""
        header["UTCTIME"] = header["UTCISO"].replace(" ", "T")
        img[0].header = header
        return img[0].data, img[0].header

    # def apply_reduction(self, raw_image_list):
    #     return
