from winterdrp.processors.database.database_exporter import DatabaseImageExporter
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp, PSFex
from winterdrp.pipelines.summer.config import get_summer_schema_path, summer_weight_path, \
    sextractor_astrometry_config, sextractor_photometry_config, scamp_path, swarp_config_path, psfex_config_path, \
    sextractor_candidates_config, PIPELINE_NAME, SUMMER_PIXEL_SCALE, summer_cal_requirements, summer_mask_path
from winterdrp.pipelines.summer.config.schema import summer_schema_dir
from winterdrp.processors.utils import ImageSaver, ImageLoader, ImageSelector, ImageBatcher
from winterdrp.processors.utils.cal_hunter import CalHunter
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
from winterdrp.processors.mask import MaskPixels
from winterdrp.pipelines.summer.load_summer_image import load_raw_summer_image, load_proc_summer_image
from winterdrp.pipelines.summer.generator import summer_astrometric_catalog_generator, \
    summer_photometric_catalog_generator, summer_reference_image_generator, summer_reference_psfex, \
    summer_reference_image_resampler, summer_reference_sextractor


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
        db_name=PIPELINE_NAME,
        db_table="exposures",
        schema_path=get_summer_schema_path("exposures"),
        full_setup=True,
        schema_dir=summer_schema_dir
    ),
    MaskPixels(mask_path=summer_mask_path),
    DatabaseImageExporter(
        db_name=PIPELINE_NAME,
        db_table="raw",
        schema_path=get_summer_schema_path("raw"),
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
    AutoAstrometry(pa=0, inv=True, pixel_scale=SUMMER_PIXEL_SCALE),
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
        db_name=PIPELINE_NAME,
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
        db_name=PIPELINE_NAME,
        db_table="diff",
        schema_path=get_summer_schema_path("diff"),
    ),
]

extract_candidates = [
    DetectCandidates(
        output_sub_dir="subtract",
        **sextractor_candidates_config
    ),
    RegionsWriter(output_dir_name='candidates'),
    PSFPhotometry(),
    AperturePhotometry(aper_diameters=[8, 40], cutout_size_aper_phot=100, bkg_in_diameters=[25, 90],
                       bkg_out_diameters=[40, 100], col_suffix_list=['', 'big']),
    DataframeWriter(output_dir_name='candidates'),
]

imsub = subtract + export_diff_to_db + export_diff_to_db