import os
import astropy.io.fits
import numpy as np
from astropy.io.fits import HDUList
from winterdrp.pipelines.base_pipeline import Pipeline

from winterdrp.processors.dark import DarkCalibrator
from winterdrp.processors.flat import SkyFlatCalibrator, OldSkyFlatCalibrator
from winterdrp.processors.sky import NightSkyMedianCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path, sextractor_astrometry_config, scamp_fp_path, \
    swarp_sp_path
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.catalog import Gaia2Mass


wirc_gain = 1.2


def wirc_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30)


class WircPipeline(Pipeline):

    name = "wirc"

    # photometry_cal = {
    #     "J": ()
    # }

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
            MaskPixels(mask_path=wirc_mask_path),
            DarkCalibrator(),
            SkyFlatCalibrator(),
            NightSkyMedianCalibrator(),
            ImageSaver(output_dir_name="testa"),
            AutoAstrometry(),
            ImageSaver(output_dir_name="testb"),
            # Sextractor(
            #     output_sub_dir="postprocess",
            #     **sextractor_astrometry_config
            # ),
            # ImageSaver(output_dir_name="testc"),
            # Scamp(
            #     ref_catalog_generator=wirc_astrometric_catalog_generator,
            #     scamp_config_path=scamp_fp_path,
            # ),
            # ImageSaver(output_dir_name="testd"),
            # Swarp(swarp_config_path=swarp_sp_path),
            # ImageSaver(output_dir_name="latest"),
        ]
    }

    @staticmethod
    def reformat_raw_data(
            img: HDUList,
            path: str
    ) -> [np.array, astropy.io.fits.Header]:
        header = img[0].header
        header["FILTER"] = header["AFT"].split("__")[0]
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]
        if "CALSTEPS" not in header.keys():
            header["CALSTEPS"] = ""
        header["BASENAME"] = os.path.basename(path)
        header["TARGET"] = header["OBJECT"].lower()
        header["UTCTIME"] = header["UTSHUT"]
        header.append(('GAIN', wirc_gain, 'Gain in electrons / ADU'), end=True)
        img[0].header = header
        return img[0].data, img[0].header

    # def apply_reduction(self, raw_image_list):
    #     return
