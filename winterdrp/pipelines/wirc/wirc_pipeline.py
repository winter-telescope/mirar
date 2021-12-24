import os

import astropy.io.fits
import numpy as np
from astropy.io.fits import HDUList
from winterdrp.pipelines.base_pipeline import Pipeline

wirc_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class WircPipeline(Pipeline):

    name = "wirc"

    astrometry_cal = ("GAIA", 9., 13.)
    photometry_cal = {
        "J": ()
    }

    # Set up elements to use

    header_keys = [
        "UTSHUT",
        'OBJECT',
        "FILTER",
        "EXPTIME",
        "COADDS",
    ]

    batch_split_keys = ["OBJECT", "FILTER"]

    pipeline_configurations = {
        None: [
            ("dark",),
            # "flat",
            ("save", "preprocess"),
            # "sextractor",
            # "stack",
            # "dither"
        ]
    }

    @staticmethod
    def reformat_raw_data(
            img: HDUList
    ) -> [np.array, astropy.io.fits.Header]:
        header = img[0].header
        header["FILTER"] = header["AFT"].split("__")[0]
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]
        header["CALSTEPS"] = ""
        img[0].header = header
        return img[0].data, img[0].header

    # def apply_reduction(self, raw_image_list):
    #     return
