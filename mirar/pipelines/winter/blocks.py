"""
Module for WINTER data reduction
"""
# pylint: disable=duplicate-code
from mirar.catalog.kowalski import PS1, TMASS
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import (
    BASE_NAME_KEY,
    CAND_NAME_KEY,
    DITHER_N_KEY,
    EXPTIME_KEY,
    MAX_DITHER_KEY,
    OBSCLASS_KEY,
    TARGET_KEY,
    ZP_KEY,
    base_output_dir,
)
from mirar.pipelines.winter.config import (
    prv_candidate_cols,
    psfex_path,
    scamp_config_path,
    sextractor_anet_config,
    sextractor_astrometry_config,
    sextractor_astromstats_config,
    sextractor_candidate_config,
    sextractor_photometry_config,
    sextractor_reference_config,
    swarp_config_path,
    winter_avro_schema_path,
)
from mirar.pipelines.winter.constants import NXSPLIT, NYSPLIT
from mirar.pipelines.winter.generator import (
    winter_astrometric_ref_catalog_generator,
    winter_astrometric_ref_catalog_namer,
    winter_astrometry_sextractor_catalog_purifier,
    winter_astrostat_catalog_purifier,
    winter_candidate_annotator_filterer,
    winter_candidate_avro_fields_calculator,
    winter_fourier_filtered_image_generator,
    winter_history_deprecated_constraint,
    winter_photometric_catalog_generator,
    winter_photometric_ref_catalog_namer,
    winter_reference_generator,
    winter_reference_image_resampler_for_zogy,
    winter_reference_psfex,
    winter_reference_sextractor,
    winter_stackid_annotator,
)
from mirar.pipelines.winter.load_winter_image import (
    annotate_winter_subdet_headers,
    get_raw_winter_mask,
    load_stacked_winter_image,
    load_test_stacked_winter_image,
    load_test_winter_image,
    load_winter_mef_image,
)
from mirar.pipelines.winter.models import (
    CANDIDATE_PREFIX,
    DEFAULT_FIELD,
    NAME_START,
    AstrometryStat,
    Candidate,
    Diff,
    Exposure,
    Raw,
    Stack,
)
from mirar.pipelines.winter.validator import (
    masked_images_rejector,
    poor_astrometric_quality_rejector,
)
from mirar.processors.astromatic import PSFex, Scamp
from mirar.processors.astromatic.sextractor.background_subtractor import (
    SextractorBkgSubtractor,
)
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.astrometry.validate import AstrometryStatsWriter
from mirar.processors.avro import IPACAvroExporter
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_inserter import (
    DatabaseImageBatchInserter,
    DatabaseImageInserter,
    DatabaseSourceInserter,
)
from mirar.processors.database.database_selector import (
    CrossmatchSourceWithDatabase,
    DatabaseHistorySelector,
)
from mirar.processors.database.database_updater import ImageDatabaseMultiEntryUpdater
from mirar.processors.mask import (  # MaskAboveThreshold,
    MaskDatasecPixels,
    MaskPixelsFromFunction,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.photometry.aperture_photometry import SourceAperturePhotometry
from mirar.processors.photometry.psf_photometry import SourcePSFPhotometry
from mirar.processors.reference import GetReferenceImage, ProcessReference
from mirar.processors.sky import SkyFlatCalibrator
from mirar.processors.sources import (
    CandidateNamer,
    CustomSourceModifier,
    SourceLoader,
    SourceWriter,
    ZOGYSourceDetector,
)
from mirar.processors.split import SUB_ID_KEY, SplitImage
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageSaver,
    ImageSelector,
    MEFLoader,
)
from mirar.processors.xmatch import XMatch
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

build_test = [
    MEFLoader(
        input_sub_dir="raw",
        load_image=load_winter_mef_image,
    ),
    ImageSelector(
        ("BOARD_ID", "4"),
        (OBSCLASS_KEY, ["dark", "science"]),
        (EXPTIME_KEY, "120.0"),
        ("FIELDID", ["9767", str(DEFAULT_FIELD)]),
    ),
    ImageSaver("testdata", output_dir=get_test_data_dir()),
]

load_test = [
    ImageLoader(
        input_img_dir=get_test_data_dir(),
        input_sub_dir="raw",
        load_image=load_test_winter_image,
    ),
    ImageBatcher("UTCTIME"),
]

load_ref = [
    ImageLoader(
        input_sub_dir="stack",
        load_image=load_stacked_winter_image,
        input_img_dir=base_output_dir,
    )
]

refbuild = [
    GetReferenceImage(ref_image_generator=winter_reference_generator),
    ImageSaver(output_dir_name="stacked_ref"),
]

# Start for new WINTER blocks:

# Loading

load_raw = [
    MEFLoader(
        input_sub_dir="raw",
        load_image=load_winter_mef_image,
    ),
]

extract_all = [
    ImageBatcher("UTCTIME"),
    DatabaseImageBatchInserter(db_table=Exposure, duplicate_protocol="ignore"),
    ImageSelector((OBSCLASS_KEY, ["dark", "science"])),
]

csvlog = [
    CSVLog(
        export_keys=[
            "UTCTIME",
            "PROGNAME",
            DITHER_N_KEY,
            MAX_DITHER_KEY,
            "FILTER",
            EXPTIME_KEY,
            OBSCLASS_KEY,
            "BOARD_ID",
            "BASENAME",
            TARGET_KEY,
            "RADEG",
            "DECDEG",
            "T_ROIC",
            "FIELDID",
        ]
    ),
]

select_split_subset = [ImageSelector(("SUBCOORD", "0_0"))]

# Optional subset selection
BOARD_ID = 4
select_subset = [
    ImageSelector(
        ("EXPTIME", "120.0"),
        ("FIELDID", ["3944", "999999999", "6124"]),
        ("BOARD_ID", str(BOARD_ID)),
        ("FILTER", ["dark", "J"]),
    ),
]

select_ref = [
    ImageSelector(
        ("FIELDID", str(3944)),
        ("BOARD_ID", str(BOARD_ID)),
    ),
    ImageDebatcher(),
    ImageBatcher("STACKID"),
]

# Split

mask_and_split = [
    ImageBatcher(BASE_NAME_KEY),
    # MaskAboveThreshold(threshold=40000.0),
    MaskDatasecPixels(),
    MaskPixelsFromFunction(mask_function=get_raw_winter_mask),
    SplitImage(n_x=NXSPLIT, n_y=NYSPLIT),
    CustomImageBatchModifier(annotate_winter_subdet_headers),
]

# Save raw images

save_raw = [
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
    DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace"),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", TARGET_KEY, "SUBCOORD"]),
    CustomImageBatchModifier(winter_stackid_annotator),
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Various processing steps

# Load from unpacked dir

load_unpacked = [
    ImageLoader(input_sub_dir="raw_unpacked"),
    ImageBatcher("UTCTIME"),
    CSVLog(
        export_keys=[
            "UTCTIME",
            "PROGNAME",
            DITHER_N_KEY,
            MAX_DITHER_KEY,
            "FILTER",
            EXPTIME_KEY,
            OBSCLASS_KEY,
            "BOARD_ID",
            "BASENAME",
            TARGET_KEY,
            "RADEG",
            "DECDEG",
            "T_ROIC",
            "FIELDID",
            "MEDCOUNT",
        ]
    ),
]

# Detrend blocks

dark_calibrate = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", EXPTIME_KEY, "SUBCOORD"]),
    DarkCalibrator(cache_sub_dir="calibration"),
    ImageSelector((OBSCLASS_KEY, ["science"])),
    ImageSaver(output_dir_name="darkcal"),
    ImageDebatcher(),
]

flat_calibrate = [
    ImageBatcher(["BOARD_ID", "FILTER", EXPTIME_KEY, "SUBCOORD"]),
    SkyFlatCalibrator(cache_sub_dir="skycals"),
    # ImageSaver(output_dir_name="skyflatcal"),
    ImageBatcher(["BOARD_ID", "UTCTIME", "SUBCOORD"]),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    ImageSaver(output_dir_name="skysub"),
]

fourier_filter = [CustomImageBatchModifier(winter_fourier_filtered_image_generator)]

astrometry = [
    ImageDebatcher(),
    ImageBatcher(["UTCTIME", "BOARD_ID", "SUBCOORD"]),
    # ImageSaver(output_dir_name="pre_anet"),
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[15, 23],
        scale_units="amw",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=1.0,
        sextractor_config_path=sextractor_anet_config["config_path"],
        use_weight=True,
        timeout=120,
    ),
    ImageSaver(output_dir_name="post_anet"),
    ImageDebatcher(),
    ImageBatcher([TARGET_KEY, "FILTER", EXPTIME_KEY, "BOARD_ID", "SUBCOORD"]),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="scamp",
        catalog_purifier=winter_astrometry_sextractor_catalog_purifier,
    ),
    CustomImageBatchModifier(winter_astrometric_ref_catalog_namer),
    Scamp(
        scamp_config_path=scamp_config_path,
        ref_catalog_generator=winter_astrometric_ref_catalog_generator,
        copy_scamp_header_to_image=True,
    ),
    ImageSaver(output_dir_name="post-scamp"),
]

validate_astrometry = [
    ImageDebatcher(),
    ImageBatcher(["UTCTIME", "BOARD_ID", "SUBCOORD"]),
    Sextractor(
        **sextractor_astromstats_config,
        write_regions_bool=True,
        output_sub_dir="astrostats",
    ),
    AstrometryStatsWriter(
        ref_catalog_generator=winter_astrometric_ref_catalog_generator,
        image_catalog_purifier=winter_astrostat_catalog_purifier,
        write_regions=True,
        cache=True,
        crossmatch_radius_arcsec=5.0,
    ),
    DatabaseImageInserter(db_table=AstrometryStat, duplicate_protocol="ignore"),
    CustomImageBatchModifier(poor_astrometric_quality_rejector),
]

stack_dithers = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", TARGET_KEY, "SUBCOORD"]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=False,
        cache=False,
        center_type="ALL",
        temp_output_sub_dir="stack_all",
        header_keys_to_combine=["RAWID"],
    ),
    ImageSaver(output_dir_name="stack"),
    CustomImageBatchModifier(masked_images_rejector),
]

photcal_and_export = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
    ),
    CustomImageBatchModifier(winter_photometric_ref_catalog_namer),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        temp_output_sub_dir="phot",
        write_regions=True,
        cache=True,
    ),
    ImageSaver(output_dir_name="final"),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="replace"),
    ImageDatabaseMultiEntryUpdater(
        sequence_key="rawid",
        db_table=Raw,
        db_alter_columns="ustackid",
    ),
]

# Image subtraction

load_test_stack = [
    ImageLoader(
        input_img_dir=get_test_data_dir(),
        input_sub_dir="final",
        load_image=load_test_stacked_winter_image,
    ),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="replace"),
]

load_stack = [
    ImageLoader(input_sub_dir="final", input_img_dir=base_output_dir),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
    # ImageSelector(
    #     (BASE_NAME_KEY, "WINTERcamera_20230727-035357-778_mef_4_0_1_stack.fits")
    # ),
]

imsub = [
    HeaderAnnotator(input_keys=[SUB_ID_KEY], output_key="SUBDETID"),
    ProcessReference(
        ref_image_generator=winter_reference_generator,
        swarp_resampler=winter_reference_image_resampler_for_zogy,
        sextractor=winter_reference_sextractor,
        ref_psfex=winter_reference_psfex,
    ),
    Sextractor(**sextractor_reference_config, output_sub_dir="subtract", cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="subtract", norm_fits=True),
    # ImageSaver(output_dir_name="presubtract"),
    ZOGYPrepare(
        output_sub_dir="subtract", sci_zp_header_key="ZP_AUTO", ref_zp_header_key=ZP_KEY
    ),
    # ImageSaver(output_dir_name="prezogy"),
    ZOGY(
        output_sub_dir="subtract", sci_zp_header_key="ZP_AUTO", ref_zp_header_key=ZP_KEY
    ),
    ImageSaver(output_dir_name="diffs"),
    DatabaseImageInserter(db_table=Diff, duplicate_protocol="replace"),
]

detect_candidates = [
    ZOGYSourceDetector(
        output_sub_dir="subtract",
        **sextractor_candidate_config,
        copy_image_keywords=["stackid", "progname"],
    ),
    SourcePSFPhotometry(),
    SourceAperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    CustomSourceModifier(winter_candidate_annotator_filterer),
    SourceWriter(output_dir_name="candidates"),
]
#
# candidate_colnames = get_column_names_from_schema(winter_candidate_config)

load_candidates = [
    SourceLoader(input_dir_name="candidates"),
]

process_candidates = [
    # SourceWriter(output_dir_name="candidates"),
    XMatch(catalog=TMASS(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1(num_sources=3, search_radius_arcmin=0.5)),
    SourceWriter(output_dir_name="kowalski"),
    CrossmatchSourceWithDatabase(
        db_table=Candidate,
        db_output_columns=[CAND_NAME_KEY],
        crossmatch_radius_arcsec=2.0,
        max_num_results=1,
    ),
    CandidateNamer(
        db_table=Candidate,
        base_name=CANDIDATE_PREFIX,
        name_start=NAME_START,
    ),
    DatabaseHistorySelector(
        crossmatch_radius_arcsec=2.0,
        time_field_name="jd",
        history_duration_days=500.0,
        db_table=Candidate,
        db_output_columns=prv_candidate_cols + [CAND_NAME_KEY],
        additional_query_constraints=winter_history_deprecated_constraint,
    ),
    CustomSourceModifier(modifier_function=winter_candidate_avro_fields_calculator),
    DatabaseSourceInserter(
        db_table=Candidate,
        duplicate_protocol="fail",
    ),
    SourceWriter(output_dir_name="preavro"),
    IPACAvroExporter(
        topic_prefix="winter",
        base_name="WNTR",
        broadcast=False,
        avro_schema_path=winter_avro_schema_path,
    ),
]

# Combinations of different blocks, to be used in configurations
process_and_stack = astrometry + validate_astrometry + stack_dithers

unpack_subset = (
    load_raw + extract_all + csvlog + select_subset + mask_and_split + save_raw
)

unpack_all = load_raw + extract_all + csvlog + mask_and_split + save_raw

full_reduction = (
    dark_calibrate
    + flat_calibrate
    + fourier_filter
    + process_and_stack
    + photcal_and_export
)

photcal_stacks = [ImageLoader(input_sub_dir="stack")] + photcal_and_export

reduce_unpacked = load_unpacked + full_reduction

reduce = unpack_all + full_reduction

reftest = (
    unpack_subset
    + dark_calibrate
    + flat_calibrate
    + process_and_stack
    + select_ref
    + refbuild
)

detrend_unpacked = load_unpacked + dark_calibrate + flat_calibrate

only_ref = load_ref + select_ref + refbuild

realtime = extract_all + mask_and_split + save_raw + full_reduction

candidates = detect_candidates + process_candidates

full = realtime + imsub
