"""
Script containing the various
:class:`~mirar.processors.base_processor.BaseProcessor`
lists which are used to build configurations for the
:class:`~mirar.pipelines.summer.summer_pipeline.SummerPipeline`.
"""

# pylint: disable=duplicate-code
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import BASE_NAME_KEY, GAIN_KEY, OBSCLASS_KEY, core_fields
from mirar.pipelines.summer.config import (
    SUMMER_PIXEL_SCALE,
    psfex_config_path,
    scamp_path,
    sextractor_astrometry_config,
    sextractor_candidates_config,
    sextractor_photometry_config,
    summer_cal_requirements,
    summer_mask_path,
    swarp_config_path,
)
from mirar.pipelines.summer.generator import (
    summer_astrometric_catalog_generator,
    summer_photometric_catalog_generator,
    summer_reference_image_generator,
    summer_reference_image_resampler,
    summer_reference_psfex,
    summer_reference_sextractor,
    summer_zogy_catalogs_purifier,
)
from mirar.pipelines.summer.load_summer_image import (
    load_proc_summer_image,
    load_raw_summer_image,
)
from mirar.pipelines.summer.models import Diff, Exposure, Proc, Raw
from mirar.processors import BiasCalibrator, FlatCalibrator
from mirar.processors.astromatic import PSFex, Scamp, Sextractor, Swarp
from mirar.processors.astrometry.autoastrometry import AutoAstrometry
from mirar.processors.cosmic_rays import LACosmicCleaner
from mirar.processors.csvlog import CSVLog
from mirar.processors.database.database_inserter import DatabaseImageInserter
from mirar.processors.database.database_updater import ImageSequenceDatabaseUpdater
from mirar.processors.mask import MaskPixelsFromPath
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import ProcessReference
from mirar.processors.sources import SourceWriter
from mirar.processors.sources.source_detector import ZOGYSourceDetector
from mirar.processors.sources.utils import RegionsWriter
from mirar.processors.utils import ImageBatcher, ImageLoader, ImageSaver, ImageSelector
from mirar.processors.utils.cal_hunter import CalHunter
from mirar.processors.utils.header_annotate import HeaderEditor
from mirar.processors.utils.simulate_realtime import RealtimeImageSimulator
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [
    ImageLoader(load_image=load_raw_summer_image),
]

load_test = [
    ImageLoader(
        input_img_dir=get_test_data_dir(),
        input_sub_dir="raw",
        load_image=load_raw_summer_image,
    ),
]

load_test_proc = [
    ImageLoader(
        input_img_dir=get_test_data_dir(),
        input_sub_dir="processed",
        load_image=load_proc_summer_image,
    ),
    ImageSelector((BASE_NAME_KEY, "SUMMER_20220816_042349_Camera0.resamp.fits")),
]

sim_realtime = [
    RealtimeImageSimulator(
        input_img_dir=get_test_data_dir(),
        input_img_names=[
            "summer/20220402/raw/SUMMER_20220402_193104_Camera0.fits",
            "summer/20220402/raw/SUMMER_20220402_214324_Camera0.fits",
        ],
        output_dir=get_test_data_dir(),
        output_dir_name="raw",
    )
]

build_log = [
    CSVLog(
        export_keys=[
            "UTC",
            "FIELDID",
            "FILTERID",
            "EXPTIME",
            OBSCLASS_KEY,
            "RA",
            "DEC",
            "TARGTYPE",
            "PROGID",
            "PROGPI",
            BASE_NAME_KEY,
        ]
        + core_fields
    ),
]

load_processed = [
    ImageLoader(input_sub_dir="processed", load_image=load_proc_summer_image)
]

export_raw = [
    # ImageSelector(("BASENAME", "SUMMER_20220402_214324_Camera0.fits")),
    DatabaseImageInserter(
        db_table=Exposure,
        duplicate_protocol="replace",
    ),
    MaskPixelsFromPath(mask_path=summer_mask_path),
    DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace"),
    ImageSelector((OBSCLASS_KEY, ["bias", "flat", "science"])),
]

cal_hunter = [
    CalHunter(load_image=load_raw_summer_image, requirements=summer_cal_requirements),
]

test_cr = [
    MaskPixelsFromPath(mask_path=summer_mask_path),
    BiasCalibrator(),
    ImageSelector((OBSCLASS_KEY, ["bias", "science"])),
    ImageBatcher(split_key="filter"),
    FlatCalibrator(),
    ImageSelector((OBSCLASS_KEY, ["science"])),
    LACosmicCleaner(effective_gain_key=GAIN_KEY, readnoise=2),
    ImageSaver(output_dir_name="crclean"),
]

process_raw = [
    BiasCalibrator(),
    ImageSelector((OBSCLASS_KEY, ["flat", "science"])),
    ImageBatcher(split_key="filter"),
    FlatCalibrator(),
    ImageBatcher(split_key=BASE_NAME_KEY),
    ImageSelector((OBSCLASS_KEY, ["science"])),
    LACosmicCleaner(effective_gain_key=GAIN_KEY, readnoise=2),
    ImageSaver(output_dir_name="detrend", write_mask=True),
    AutoAstrometry(pa=0, inv=True, pixel_scale=SUMMER_PIXEL_SCALE),
    ImageSaver(output_dir_name="detrend", write_mask=True),
    Sextractor(
        output_sub_dir="sextractor",
        # TODO: work out why this was ever here...
        # weight_image=summer_weight_path,
        checkimage_name=None,
        checkimage_type=None,
        **sextractor_astrometry_config,
    ),
    Scamp(
        ref_catalog_generator=summer_astrometric_catalog_generator,
        scamp_config_path=scamp_path,
        cache=True,
    ),
    Swarp(
        swarp_config_path=swarp_config_path,
        cache=True,
        # TODO: work out why this was ever here...
        # imgpixsize=2400
    ),
    ImageSaver(output_dir_name="processed", write_mask=True),
    Sextractor(
        output_sub_dir="photprocess",
        checkimage_type="BACKGROUND_RMS",
        **sextractor_photometry_config,
    ),
    PhotCalibrator(
        ref_catalog_generator=summer_photometric_catalog_generator,
    ),
    ImageSaver(
        output_dir_name="processed",
        # TODO: work out why this was ever here...
        # additional_headers=["PROCIMG"],
        write_mask=True,
    ),
    HeaderEditor(edit_keys="procstatus", values=1),
    DatabaseImageInserter(db_table=Proc, duplicate_protocol="replace"),
    ImageSequenceDatabaseUpdater(
        # db_name=DB_NAME,
        db_table=Raw,
        db_alter_columns="procstatus",
    ),
]

# standard_summer_reduction = export_raw + cal_hunter + process_raw #FIXME


subtract = [
    ImageBatcher(split_key=BASE_NAME_KEY),
    ImageSelector((OBSCLASS_KEY, "science")),
    ProcessReference(
        ref_image_generator=summer_reference_image_generator,
        ref_psfex=summer_reference_psfex,
        sextractor=summer_reference_sextractor,
        swarp_resampler=summer_reference_image_resampler,
    ),
    Sextractor(
        output_sub_dir="subtract",
        cache=False,
        write_regions_bool=True,
        **sextractor_photometry_config,
    ),
    PSFex(config_path=psfex_config_path, output_sub_dir="subtract", norm_fits=True),
    ImageSaver(output_dir_name="ref"),
    ZOGYPrepare(
        output_sub_dir="subtract",
        sci_zp_header_key="ZP_AUTO",
        ref_zp_header_key="ZP",
        catalog_purifier=summer_zogy_catalogs_purifier,
    ),
    ZOGY(output_sub_dir="subtract"),
]

export_diff_to_db = [
    DatabaseImageInserter(
        db_table=Diff,
    ),
]

extract_candidates = [
    ZOGYSourceDetector(output_sub_dir="subtract", **sextractor_candidates_config),
    RegionsWriter(output_dir_name="candidates"),
    PSFPhotometry(),
    AperturePhotometry(
        aper_diameters=[8, 40],
        phot_cutout_half_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    SourceWriter(output_dir_name="candidates"),
]

imsub = subtract + export_diff_to_db + extract_candidates  # FIXME
