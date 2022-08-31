import logging
import os
import astropy.io.fits
import numpy as np

from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.downloader.caltech import download_via_ssh
from winterdrp.processors.database.database_exporter import DatabaseImageExporter
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp, PSFex
from winterdrp.pipelines.summer.summer_files import get_summer_schema_path, summer_weight_path, \
    sextractor_astrometry_config, sextractor_photometry_config, scamp_path, swarp_config_path
from winterdrp.pipelines.summer.summer_files.schema import summer_schema_dir
from winterdrp.processors.utils import ImageSaver, ImageLoader, ImageSelector, ImageBatcher
from winterdrp.processors.utils.cal_hunter import CalHunter, CalRequirement
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors import BiasCalibrator, FlatCalibrator
from winterdrp.processors.csvlog import CSVLog
from winterdrp.paths import core_fields, base_name_key
from winterdrp.processors.reference import Reference
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare, default_summer_catalog_purifier
from winterdrp.processors.candidates.candidate_detector import DetectCandidates
from winterdrp.processors.photometry.psf_photometry import PSFPhotometry
from winterdrp.processors.photometry.aperture_photometry import AperturePhotometry
from winterdrp.processors.candidates.utils import RegionsWriter, DataframeWriter
from winterdrp.pipelines.summer.load_summer_image import load_raw_summer_image, load_proc_summer_image
from winterdrp.pipelines.summer.generator import summer_astrometric_catalog_generator, \
    summer_photometric_catalog_generator, summer_reference_image_generator, summer_reference_psfex, \
    summer_reference_image_resampler, summer_reference_sextractor
from winterdrp.pipelines.summer.summer_files import psfex_config_path

summer_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
summer_gain = 1.0
summer_pixel_scale = 0.466

logger = logging.getLogger(__name__)

summer_cal_requirements = [
    CalRequirement(target_name="bias", required_field="EXPTIME", required_values=["0.0"]),
    CalRequirement(target_name="flat", required_field="FILTERID", required_values=["u", "g", "r", "i"]),
]


pipeline_name = "summer"


load_raw = [ImageLoader(load_image=load_raw_summer_image)]

load_processed = [ImageLoader(input_sub_dir='processed', load_image=load_proc_summer_image)]

standard_summer_reduction = [
    CSVLog(
        export_keys=[
            "UTC", 'FIELDID', "FILTERID", "EXPTIME", "OBSTYPE", "RA", "DEC", "TARGTYPE","PROGID", "PROGPI",
            base_name_key
        ] + core_fields
    ),
    DatabaseImageExporter(
        db_name=pipeline_name,
        db_table="exposures",
        schema_path=get_summer_schema_path("exposures"),
        full_setup=True,
        schema_dir=summer_schema_dir
    ),
    ImageSelector(("OBSTYPE", ["BIAS", "FLAT", "SCIENCE"])),
    CalHunter(
        load_image=load_raw_summer_image,
        requirements=summer_cal_requirements
    ),
    BiasCalibrator(),
    ImageSelector(("OBSTYPE", ["FLAT", "SCIENCE"])),
    ImageBatcher(split_key="filter"),
    FlatCalibrator(),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageBatcher(base_name_key),
    AutoAstrometry(pa=0, inv=True, pixel_scale=summer_pixel_scale),
    Sextractor(
        output_sub_dir="sextractor",
        weight_image=summer_weight_path,
        checkimage_name=None,
        checkimage_type=None,
        **sextractor_astrometry_config
    ),
    Scamp(
        ref_catalog_generator=summer_astrometric_catalog_generator,
        scamp_config_path=scamp_path,
    ),
    Swarp(swarp_config_path=swarp_config_path, imgpixsize=2400),
    Sextractor(output_sub_dir="photprocess",
               checkimage_type='BACKGROUND_RMS',
               **sextractor_photometry_config),
    PhotCalibrator(ref_catalog_generator=summer_photometric_catalog_generator),
    ImageSaver(output_dir_name="processed", additional_headers=['PROCIMG'], write_mask=True),
    DatabaseImageExporter(
        db_name=pipeline_name,
        db_table="proc",
        schema_path=get_summer_schema_path("proc")
    )
]

subtract = [
    ImageBatcher(split_key=base_name_key),
    ImageSelector(('OBSTYPE', 'SCIENCE')),
    # ImageSelector((base_name_key, ["SUMMER_20220816_042926_Camera0.resamp.fits","SUMMER_20220816_042743_Camera0.resamp.fits"])),
    Reference(
        ref_image_generator=summer_reference_image_generator,
        ref_psfex=summer_reference_psfex,
        ref_sextractor=summer_reference_sextractor,
        ref_swarp_resampler=summer_reference_image_resampler
    ),
    Sextractor(
        output_sub_dir='subtract',
        cache=False,
        write_regions_file=True,
        **sextractor_photometry_config
    ),
    PSFex(config_path=psfex_config_path,
          output_sub_dir="subtract",
          norm_fits=True),
    ImageSaver(output_dir_name='ref'),
    ZOGYPrepare(output_sub_dir="subtract", sci_zp_header_key='ZP_AUTO',
                catalog_purifier=default_summer_catalog_purifier),
    ZOGY(output_sub_dir="subtract"),
]

export_diff_to_db = [
    DatabaseImageExporter(
        db_name=pipeline_name,
        db_table="diff",
        schema_path=get_summer_schema_path("diff"),
    ),
]

extract_candidates = [
    DetectCandidates(
        output_sub_dir="subtract",
        cand_det_sextractor_config='winterdrp/pipelines/summer/summer_imsub_files/config/photomCat.sex',  # Delete
        cand_det_sextractor_nnw='winterdrp/pipelines/summer/summer_imsub_files/config/default.nnw',
        cand_det_sextractor_filter='winterdrp/pipelines/summer/summer_imsub_files/config/default.conv',
        cand_det_sextractor_params='winterdrp/pipelines/summer/summer_imsub_files/config/Scorr.param'
    ),
    RegionsWriter(output_dir_name='candidates'),
    PSFPhotometry(),
    AperturePhotometry(aper_diameters=[8, 40], cutout_size_aper_phot=100, bkg_in_diameters=[25, 90],
                       bkg_out_diameters=[40, 100], col_suffix_list=['', 'big']),
    DataframeWriter(output_dir_name='candidates'),
]

imsub = subtract + export_diff_to_db + export_diff_to_db


class SummerPipeline(Pipeline):

    name = pipeline_name
    default_cal_requirements = summer_cal_requirements

    all_pipeline_configurations = {
        None: load_raw + standard_summer_reduction,
        'imsub': load_processed + imsub,
        "full": load_raw + standard_summer_reduction + imsub,
        "realtime": standard_summer_reduction,
    }

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        download_via_ssh(
            server="jagati.caltech.edu",
            base_dir="/data/viraj/winter_data/commissioning/raw/",
            night=night,
            pipeline=pipeline_name
        )

    @staticmethod
    def load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        return load_raw_summer_image(path)
