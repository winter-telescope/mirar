import os
import astropy.io.fits
import numpy as np
from astropy.io import fits
from astropy.time import Time
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.processors.dark import DarkCalibrator
from winterdrp.processors.flat import SkyFlatCalibrator
from winterdrp.processors.sky import NightSkyMedianCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path, sextractor_astrometry_config, scamp_fp_path, \
    swarp_sp_path
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.catalog import Gaia2Mass
from winterdrp.downloader.caltech import download_via_ssh
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
from winterdrp.paths import coadd_key, proc_history_key
from winterdrp.processors.csvlog import CSVLog
import logging

logger = logging.getLogger(__name__)


def wirc_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30)


def wirc_photometric_catalog_generator(
        header: astropy.io.fits.Header
):
    filter_name = header['FILTER']
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30, filter_name=filter_name)


pipeline_name = "wirc"


def load_raw_wirc_image(
        path: str
) -> tuple[np.array, astropy.io.fits.Header]:
    with fits.open(path) as img:
        data = img[0].data
        header = img[0].header
        header["FILTER"] = header["AFT"].split("__")[0]

        if header["OBJECT"] in ["acquisition", "pointing", "focus", "none"]:
            header["OBSTYPE"] = "calibration"

        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]

        header["CALSTEPS"] = ""

        header["BASENAME"] = os.path.basename(path)
        header["TARGET"] = header["OBJECT"].lower()
        header["UTCTIME"] = header["UTSHUT"]
        header["MJD-OBS"] = Time(header['UTSHUT']).mjd
        if coadd_key not in header.keys():
            logger.debug(f"No {coadd_key} entry. Setting coadds to 1.")
            header[coadd_key] = 1
        if proc_history_key not in header.keys():
            header[proc_history_key] = ""

        filter_dict = {'J': 1,'H': 2, 'Ks': 3}
        if "FILTERID" not in header.keys():
            header["FILTERID"] = filter_dict[header["FILTER"]]
        if "FIELDID" not in header.keys():
            header["FIELDID"] = 99999
        if "PROGPI" not in header.keys():
            header["PROGPI"] = "Kasliwal"
        if "PROGID" not in header.keys():
            header["PROGID"] = 0

        data = data.astype(float)
        data[data == 0.] = np.nan
    return data, header


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
            ImageLoader(
                input_sub_dir="raw",
                load_image=load_raw_wirc_image
            ),
            CSVLog(
                export_keys=["OBJECT", "FILTER", "UTSHUT", "EXPTIME", "COADDS", "OBSTYPE", "OBSCLASS"]
            ),
            MaskPixels(mask_path=wirc_mask_path),
            ImageBatcher(split_key="exptime"),
            DarkCalibrator(),
            ImageDebatcher(),
            ImageSelector(("obsclass", "science")),
            ImageBatcher(split_key="filter"),
            SkyFlatCalibrator(),
            NightSkyMedianCalibrator(),
            ImageSelector(
                ("object", "ZTF21aagppzg"),
                ("filter", "J")
            ),
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
            Sextractor(
                output_sub_dir="final_sextractor",
                **sextractor_astrometry_config
            ),
            ImageSaver(output_dir_name="final"),
            PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator),
            # ImageDebatcher(),
            # ImageBatcher(split_key="filter"),
            # SkyFlatCalibrator(),
            # NightSkyMedianCalibrator(),
            # ImageSelector(
            #     ("object", "ZTF21acbnfos"),
            # ),
            # AutoAstrometry(catalog="tmc"),
            # Sextractor(
            #     output_sub_dir="postprocess",
            #     **sextractor_astrometry_config
            # ),
            # Scamp(
            #     ref_catalog_generator=wirc_astrometric_catalog_generator,
            #     scamp_config_path=scamp_fp_path,
            # ),
            # Swarp(swarp_config_path=swarp_sp_path),
            # Sextractor(
            #     output_sub_dir="final_sextractor",
            #     **sextractor_astrometry_config
            # ),
            # ImageSaver(output_dir_name="final"),
            # PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator),
        ]
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
