import os
import astropy.io.fits
import numpy as np
from astropy.io import fits
from astropy.time import Time
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
from winterdrp.downloader.caltech import download_via_ssh


def wirc_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30)


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

    batch_split_keys = ["OBJECT", "FILTER"]

    pipeline_configurations = {
        None: [
            MaskPixels(mask_path=wirc_mask_path),
            DarkCalibrator(),
            SkyFlatCalibrator(),
            NightSkyMedianCalibrator(),
            AutoAstrometry(catalog="tmc"),
            Sextractor(
                output_sub_dir="postprocess",
                **sextractor_astrometry_config
            ),
            Scamp(
                ref_catalog_generator=wirc_astrometric_catalog_generator,
                scamp_config_path=scamp_fp_path,
            ),
            Swarp(swarp_config_path=swarp_sp_path),
            ImageSaver(output_dir_name="latest"),
        ]
    }

    def load_raw_image(
            self,
            path: str
    ) -> tuple[np.array, astropy.io.fits.Header]:
        with fits.open(path) as img:
            data = img[0].data
            header = img[0].header
            header["FILTER"] = header["AFT"].split("__")[0]
            header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]
            header["CALSTEPS"] = ""
            header["BASENAME"] = os.path.basename(path)
            header["TARGET"] = header["OBJECT"].lower()
            header["UTCTIME"] = header["UTSHUT"]
            header["MJD-OBS"] = Time(header['UTSHUT']).mjd
            header.append(('GAIN', self.gain, 'Gain in electrons / ADU'), end=True)
            header = self.set_saturation(header)
        return data, header

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        download_via_ssh(
            server="gayatri.caltech.edu",
            base_dir="/data/sanand/WIRC/",
            night=night,
            pipeline=pipeline_name,
            server_sub_dir="raw"
        )
