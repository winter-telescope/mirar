"""
Module for WINTER data reduction
"""

# pylint: disable=duplicate-code
import os

from winterrb.model import WINTERNet

from mirar.catalog.kowalski import PS1, PS1STRM, TMASS, Gaia, GaiaBright, PS1SGSc
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.paths import (
    BASE_NAME_KEY,
    DITHER_N_KEY,
    EXPTIME_KEY,
    LATEST_SAVE_KEY,
    MAX_DITHER_KEY,
    OBSCLASS_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    SOURCE_HISTORY_KEY,
    SOURCE_NAME_KEY,
    TARGET_KEY,
    ZP_KEY,
    base_output_dir,
)
from mirar.pipelines.winter.config import (
    prv_candidate_cols,
    psfex_path,
    scamp_config_path,
    sextractor_astrometry_config,
    sextractor_astromstats_config,
    sextractor_candidate_config,
    sextractor_photometry_config,
    sextractor_photometry_psf_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
    winter_avro_schema_path,
    winter_cal_requirements,
    winter_fritz_config,
)
from mirar.pipelines.winter.constants import NXSPLIT, NYSPLIT
from mirar.pipelines.winter.generator import (
    apply_rb_to_table,
    mask_stamps_around_bright_stars,
    select_winter_sky_flat_images,
    winter_anet_sextractor_config_path_generator,
    winter_astrometric_ref_catalog_generator,
    winter_astrometric_ref_catalog_namer,
    winter_astrometry_sextractor_catalog_purifier,
    winter_astrostat_catalog_purifier,
    winter_boardid_6_demasker,
    winter_candidate_annotator_filterer,
    winter_candidate_avro_fields_calculator,
    winter_candidate_quality_filterer,
    winter_fourier_filtered_image_generator,
    winter_history_deprecated_constraint,
    winter_imsub_catalog_purifier,
    winter_new_source_updater,
    winter_photcal_color_columns_generator,
    winter_photometric_catalog_generator,
    winter_photometric_catalogs_purifier,
    winter_photometric_ref_catalog_namer,
    winter_reference_generator,
    winter_reference_image_resampler_for_zogy,
    winter_reference_psf_phot_sextractor,
    winter_reference_psfex,
    winter_reference_sextractor,
    winter_skyportal_annotator,
    winter_source_entry_updater,
    winter_stackid_annotator,
)
from mirar.pipelines.winter.load_winter_image import (
    annotate_winter_subdet_headers,
    get_raw_winter_mask,
    load_astrometried_winter_image,
    load_stacked_winter_image,
    load_test_winter_image,
    load_winter_mef_image,
    load_winter_stack,
)
from mirar.pipelines.winter.models import (
    DEFAULT_FIELD,
    NAME_START,
    SOURCE_PREFIX,
    AstrometryStat,
    Candidate,
    Diff,
    Exposure,
    Raw,
    Source,
    Stack,
)
from mirar.pipelines.winter.validator import (
    masked_images_rejector,
    poor_astrometric_quality_rejector,
    winter_dark_oversubtraction_rejector,
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
from mirar.processors.catalog_limiting_mag import CatalogLimitingMagnitudeCalculator
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_inserter import (
    DatabaseImageBatchInserter,
    DatabaseImageInserter,
    DatabaseSourceInserter,
)
from mirar.processors.database.database_selector import SelectSourcesWithMetadata
from mirar.processors.database.database_updater import ImageDatabaseMultiEntryUpdater
from mirar.processors.flat import FlatCalibrator
from mirar.processors.mask import (  # MaskAboveThreshold,
    MaskDatasecPixels,
    MaskPixelsFromFunction,
)
from mirar.processors.photcal import ZPWithColorTermCalculator
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import GetReferenceImage, ProcessReference
from mirar.processors.skyportal.skyportal_candidate import SkyportalCandidateUploader
from mirar.processors.sources import (
    CandidateNamer,
    CustomSourceTableModifier,
    ForcedPhotometryDetector,
    SourceBatcher,
    SourceLoader,
    SourceWriter,
    ZOGYSourceDetector,
)
from mirar.processors.sources.machine_learning import Pytorch
from mirar.processors.split import SUB_ID_KEY, SplitImage, SwarpImageSplitter
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    HeaderEditor,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImagePlotter,
    ImageRebatcher,
    ImageRejector,
    ImageSaver,
    ImageSelector,
    MEFLoader,
)
from mirar.processors.utils.cal_hunter import CalHunter
from mirar.processors.xmatch import XMatch
from mirar.processors.zogy.reference_aligner import AlignReference
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

build_test = [
    MEFLoader(
        input_sub_dir="raw",
        load_image=load_winter_mef_image,
    ),
    ImageBatcher("UTCTIME"),
    CSVLog(
        export_keys=[
            "UTCTIME",
            "PROGNAME",
            "TARGNAME",
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
            "FOCPOS",
        ]
    ),
    ImageSelector(
        ("BOARD_ID", "2"),
        (OBSCLASS_KEY, ["dark", "science"]),
        (EXPTIME_KEY, "120.0"),
        ("filter", ["dark", "J"]),
        ("FIELDID", ["3944", str(DEFAULT_FIELD)]),
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
    CalHunter(load_image=load_winter_mef_image, requirements=winter_cal_requirements),
]

load_astrometry = [
    ImageLoader(input_sub_dir="post_scamp", load_image=load_astrometried_winter_image)
]

extract_all = [
    ImageRebatcher("EXPID"),
    DatabaseImageBatchInserter(db_table=Exposure, duplicate_protocol="replace"),
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
            "MIRCOVER",
            "BASENAME",
            TARGET_KEY,
            "RADEG",
            "DECDEG",
            "T_ROIC",
            "FIELDID",
            "READOUTM",
            "READOUTV",
        ]
    ),
    ImageRebatcher(BASE_NAME_KEY),
]

select_split_subset = [ImageSelector(("SUBCOORD", "0_0"))]

# Optional subset selection
BOARD_ID = 4
select_subset = [
    ImageSelector(
        ("BOARD_ID", str(BOARD_ID)),
    ),
]

select_ref = [
    ImageSelector(
        ("FIELDID", str(3944)),
        ("BOARD_ID", str(BOARD_ID)),
    ),
    ImageRebatcher("STACKID"),
]

# mask
mask = [
    ImageSelector((OBSCLASS_KEY, ["dark", "science", "flat"])),
    # MaskAboveThreshold(threshold=40000.0),
    MaskDatasecPixels(),
    MaskPixelsFromFunction(mask_function=get_raw_winter_mask),
]

# Split
split = [
    SplitImage(n_x=NXSPLIT, n_y=NYSPLIT),
    CustomImageBatchModifier(annotate_winter_subdet_headers),
]

mask_and_split = mask + split

# Save raw images

save_raw = [
    # Group into planned stacks, and label each image with the intended stackid
    ImageRebatcher(["BOARD_ID", "FILTER", "EXPTIME", TARGET_KEY, "SUBCOORD"]),
    CustomImageBatchModifier(winter_stackid_annotator),
    # Process each raw image in parallel
    ImageRebatcher(BASE_NAME_KEY),
    ImageRejector(("BOARD_ID", "0")),
    # HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
    DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace"),
]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Various processing steps

# Load from unpacked dir

load_unpacked = [
    ImageLoader(input_sub_dir="raw_unpacked", input_img_dir=base_output_dir),
    ImageRebatcher("EXPID"),
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
    ImageRebatcher(BASE_NAME_KEY),
]

export_unpacked = [DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace")]
load_and_export_unpacked = load_unpacked + export_unpacked


# Detrend blocks

dark_calibrate = [
    ImageRebatcher(
        ["BOARD_ID", EXPTIME_KEY, "SUBCOORD", "GAINCOLT", "GAINCOLB", "GAINROW"]
    ),
    DarkCalibrator(
        cache_sub_dir="calibration_darks",
        cache_image_name_header_keys=[EXPTIME_KEY, "BOARD_ID"],
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="darkcal"),
    ImageSelector((OBSCLASS_KEY, ["science", "flat"])),
    CustomImageBatchModifier(winter_dark_oversubtraction_rejector),
]

flat_calibrate = [
    ImageSelector((OBSCLASS_KEY, ["science"])),
    ImageRebatcher(
        [
            "BOARD_ID",
            "FILTER",
            "SUBCOORD",
            "GAINCOLT",
            "GAINCOLB",
            "GAINROW",
            TARGET_KEY,
        ]
    ),
    FlatCalibrator(
        cache_sub_dir="sky_dither_flats",
        select_flat_images=select_winter_sky_flat_images,
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="skyflatcal"),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    ImageSaver(output_dir_name="skysub"),
]

load_calibrated = [
    ImageLoader(input_sub_dir="skysub", input_img_dir=base_output_dir),
    ImageBatcher(["UTCTIME", "BOARD_ID"]),
]

fourier_filter = [
    CustomImageBatchModifier(winter_fourier_filtered_image_generator),
]

astrometry = [
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[1.0, 1.3],
        scale_units="app",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=5.0,
        sextractor_config_path=winter_anet_sextractor_config_path_generator,
        use_weight=True,
        timeout=120,
        cache=False,
    ),
    ImageSaver(output_dir_name="post_anet"),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="scamp",
        catalog_purifier=winter_astrometry_sextractor_catalog_purifier,
    ),
    CustomImageBatchModifier(winter_astrometric_ref_catalog_namer),
    ImageRebatcher([TARGET_KEY, "FILTER", EXPTIME_KEY, "BOARD_ID", "SUBCOORD"]),
    Scamp(
        scamp_config_path=scamp_config_path,
        ref_catalog_generator=winter_astrometric_ref_catalog_generator,
        copy_scamp_header_to_image=True,
        cache=True,
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="post_scamp"),
]

validate_astrometry = [
    Sextractor(
        **sextractor_astromstats_config,
        write_regions_bool=True,
        output_sub_dir="astrostats",
    ),
    AstrometryStatsWriter(
        ref_catalog_generator=winter_astrometric_ref_catalog_generator,
        image_catalog_purifier=winter_astrostat_catalog_purifier,
        write_regions=True,
        cache=False,
        crossmatch_radius_arcsec=5.0,
    ),
    DatabaseImageInserter(db_table=AstrometryStat, duplicate_protocol="ignore"),
    CustomImageBatchModifier(poor_astrometric_quality_rejector),
]

stack_dithers = [
    CustomImageBatchModifier(winter_boardid_6_demasker),
    ImageRebatcher("STACKID"),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=True,
        subtract_bkg=False,
        cache=False,
        center_type="ALL",
        temp_output_sub_dir="stacks_weights",
        header_keys_to_combine=["RAWID"],
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="stack"),
]

photcal_and_export = [
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    CustomImageBatchModifier(masked_images_rejector),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="stack_psf",
        checkimage_type="BACKGROUND_RMS",
    ),
    PSFex(config_path=psfex_path, output_sub_dir="phot", norm_fits=True),
    Sextractor(
        **sextractor_photometry_psf_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
        use_psfex=True,
    ),
    CustomImageBatchModifier(winter_photometric_ref_catalog_namer),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        catalogs_purifier=winter_photometric_catalogs_purifier,
        temp_output_sub_dir="phot",
        write_regions=True,
        cache=True,
        zp_calculator=ZPWithColorTermCalculator(
            color_colnames_guess_generator=winter_photcal_color_columns_generator,
            reject_outliers=True,
            solver="curve_fit",
        ),
        zp_column_name="MAG_AUTO",
    ),
    CatalogLimitingMagnitudeCalculator(
        sextractor_mag_key_name="MAG_AUTO", write_regions=True
    ),
    AstrometryStatsWriter(
        ref_catalog_generator=winter_astrometric_ref_catalog_generator,
        image_catalog_purifier=winter_astrostat_catalog_purifier,
        write_regions=True,
        cache=False,
        crossmatch_radius_arcsec=5.0,
    ),
    ImageSaver(output_dir_name="final", compress=False),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="replace"),
    ImageDatabaseMultiEntryUpdater(
        sequence_key="rawid",
        db_table=Raw,
        db_alter_columns="ustackid",
    ),
    ImagePlotter(
        output_sub_dir="final_stacks_plots",
        annotate_fields=[
            BASE_NAME_KEY,
            "COADDS",
            TARGET_KEY,
            "CRVAL1",
            "CRVAL2",
            "FILTER",
            "ZP",
            "ZPSTD",
        ],
    ),
]

# Stack stacks together

stack_stacks = [
    ImageRebatcher(BASE_NAME_KEY),
    HeaderEditor(PROC_HISTORY_KEY, "load"),
    ImageSaver(output_dir_name="restack_masks", write_mask=True),
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
        cache=True,
    ),
    ImageSaver(output_dir_name="post_scamp"),
    HeaderAnnotator(input_keys=["TARGNAME", "FIELDID"], output_key=TARGET_KEY),
    ImageRebatcher(["SUBCOORD", "FILTER", TARGET_KEY]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=True,
        subtract_bkg=False,
        cache=False,
        center_type="ALL",
        temp_output_sub_dir="stacks_weights",
        header_keys_to_combine=["RAWID"],
    ),
    ImageSaver(output_dir_name="stack_of_stacks"),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    CustomImageBatchModifier(masked_images_rejector),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="stack_psf",
        checkimage_type="BACKGROUND_RMS",
    ),
    PSFex(config_path=psfex_path, output_sub_dir="phot", norm_fits=True),
    Sextractor(
        **sextractor_photometry_psf_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
        use_psfex=True,
    ),
    CustomImageBatchModifier(winter_photometric_ref_catalog_namer),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        catalogs_purifier=winter_photometric_catalogs_purifier,
        temp_output_sub_dir="phot",
        write_regions=True,
        cache=True,
        zp_calculator=ZPWithColorTermCalculator(
            color_colnames_guess_generator=winter_photcal_color_columns_generator,
            reject_outliers=True,
            solver="curve_fit",
        ),
        zp_column_name="MAG_AUTO",
    ),
    CatalogLimitingMagnitudeCalculator(
        sextractor_mag_key_name="MAG_AUTO", write_regions=True
    ),
    ImageSaver(output_dir_name="final_stack_of_stacks"),
]

# Image subtraction

load_final_stack = [
    ImageLoader(
        input_sub_dir="final",
        input_img_dir=base_output_dir,
        load_image=load_winter_stack,
    ),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="ignore"),
    ImageRebatcher(BASE_NAME_KEY),
]

plot_stack = [
    ImageRebatcher([TARGET_KEY, "BOARD_ID"]),
    ImagePlotter(
        output_sub_dir="final_stacks_plots",
        annotate_fields=[
            BASE_NAME_KEY,
            "COADDS",
            TARGET_KEY,
            "CRVAL1",
            "CRVAL2",
            "FILTER",
            "ZP",
            "ZPSTD",
        ],
    ),
]

split_stack = [
    ImageRebatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD", "STACKID"]),
    SwarpImageSplitter(swarp_config_path=swarp_config_path, n_x=2, n_y=1),
    ImageSaver(output_dir_name="split_stacks"),
]

imsub = [
    HeaderAnnotator(input_keys=[SUB_ID_KEY], output_key="SUBDETID"),
    ProcessReference(
        ref_image_generator=winter_reference_generator,
        swarp_resampler=winter_reference_image_resampler_for_zogy,
        sextractor=winter_reference_sextractor,
        ref_psfex=winter_reference_psfex,
        phot_sextractor=winter_reference_psf_phot_sextractor,
    ),
    Sextractor(
        **sextractor_reference_psf_phot_config,
        output_sub_dir="subtract",
        cache=False,
        use_psfex=True,
    ),
    PSFex(config_path=psfex_path, output_sub_dir="subtract", norm_fits=True),
    AlignReference(
        order=1,
        sextractor=winter_reference_sextractor,
        psfex=winter_reference_psfex,
        phot_sextractor=winter_reference_psf_phot_sextractor,
        catalog_purifier=winter_imsub_catalog_purifier,
    ),
    MaskPixelsFromFunction(mask_function=mask_stamps_around_bright_stars),
    ImageSaver(output_dir_name="presubtract"),
    ZOGYPrepare(
        output_sub_dir="subtract",
        sci_zp_header_key="ZP_AUTO",
        ref_zp_header_key=ZP_KEY,
        catalog_purifier=winter_imsub_catalog_purifier,
        x_key="XMODEL_IMAGE",
        y_key="YMODEL_IMAGE",
        flux_key="FLUX_POINTSOURCE",
    ),
    # ImageSaver(output_dir_name="prezogy"),
    ZOGY(
        output_sub_dir="subtract", sci_zp_header_key="ZP_AUTO", ref_zp_header_key=ZP_KEY
    ),
    ImageSaver(output_dir_name="diffs"),
    DatabaseImageInserter(db_table=Diff, duplicate_protocol="replace"),
    ImageSaver(output_dir_name="subtract"),
    ImageRejector(("PROGNAME", ["2023A007"])),  # Filter out galactic data
]

load_sub = [
    ImageLoader(input_sub_dir="subtract"),
]
detect_candidates = [
    ZOGYSourceDetector(
        output_sub_dir="subtract",
        **sextractor_candidate_config,
        write_regions=True,
        detect_negative_sources=True,
    ),
    PSFPhotometry(phot_cutout_half_size=10),
    AperturePhotometry(
        temp_output_sub_dir="aper_photometry",
        aper_diameters=[8, 16],
        phot_cutout_half_size=50,
        bkg_in_diameters=[25, 25],
        bkg_out_diameters=[40, 40],
        col_suffix_list=["", "big"],
    ),
    CustomSourceTableModifier(winter_candidate_annotator_filterer),
    SourceWriter(output_dir_name="candidates"),
]
#
# candidate_colnames = get_column_names_from_schema(winter_candidate_config)

load_sources = [
    SourceLoader(input_dir_name="candidates"),
    SourceBatcher(BASE_NAME_KEY),
]

ml_classify = [
    Pytorch(
        model=WINTERNet(),
        model_weights_url="https://github.com/winter-telescope/winterrb/raw/"
        "v1.0.0/models/winterrb_v1_0_0_weights.pth",
        apply_to_table=apply_rb_to_table,
    ),
    HeaderEditor(edit_keys="rbversion", values="v1.0.0"),
]

crossmatch_candidates = [
    XMatch(catalog=TMASS(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1SGSc(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1STRM(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=Gaia(num_sources=1, search_radius_arcmin=1.5)),
    XMatch(catalog=GaiaBright(num_sources=1, search_radius_arcmin=1.5)),
    CustomSourceTableModifier(
        modifier_function=winter_candidate_avro_fields_calculator
    ),
    SourceWriter(output_dir_name="kowalski"),
]

load_post_kowalski = [
    SourceLoader(input_dir_name="kowalski"),
]

select_history = [
    SelectSourcesWithMetadata(
        db_query_columns=["sourceid"],
        db_table=Candidate,
        db_output_columns=prv_candidate_cols + [SOURCE_NAME_KEY],
        base_output_column=SOURCE_HISTORY_KEY,
        additional_query_constraints=winter_history_deprecated_constraint,
    ),
]

name_candidates = [
    SourceBatcher(BASE_NAME_KEY),
    # Add the new sources to the source table
    CustomSourceTableModifier(modifier_function=winter_new_source_updater),
    # Assign names to the new sources
    CandidateNamer(
        db_table=Source,
        db_output_columns=["sourceid", SOURCE_NAME_KEY],
        base_name=SOURCE_PREFIX,
        name_start=NAME_START,
        db_name_field=SOURCE_NAME_KEY,
        crossmatch_radius_arcsec=2.0,
        ra_field_name="average_ra",
        dec_field_name="average_dec",
    ),
    # Add candidates in the candidate table
    DatabaseSourceInserter(
        db_table=Candidate,
        duplicate_protocol="fail",
    ),
    SelectSourcesWithMetadata(
        db_query_columns=["sourceid"],
        db_table=Candidate,
        db_output_columns=prv_candidate_cols + [SOURCE_NAME_KEY],
        base_output_column=SOURCE_HISTORY_KEY,
        additional_query_constraints=winter_history_deprecated_constraint,
    ),
    # Update average ra and dec for source
    CustomSourceTableModifier(modifier_function=winter_source_entry_updater),
    # Update sources in the source table
    DatabaseSourceInserter(
        db_table=Source,
        duplicate_protocol="replace",
    ),
    SourceWriter(output_dir_name="preavro"),
]

load_preavro = [
    SourceLoader(input_dir_name="preavro"),
]

avro_write = [
    # Add in the skyportal fields and all save locally
    CustomSourceTableModifier(modifier_function=winter_skyportal_annotator),
    IPACAvroExporter(
        topic_prefix="winter",
        base_name="WNTR",
        broadcast=False,
        save_local=True,
        avro_schema_path=winter_avro_schema_path,
    ),
]

# configure to broadcast to IPAC
BROADCAST_BOOL = str(os.getenv("BROADCAST_AVRO", None)) in ["True", "t", "1", "true"]

avro_broadcast = [
    # Filter out low quality candidates
    CustomSourceTableModifier(modifier_function=winter_candidate_quality_filterer),
    # Save candidates before sending to IPAC
    SourceWriter(output_dir_name="preskyportal"),
    # Only send a subset of the candidates to IPAC
    IPACAvroExporter(
        output_sub_dir="avro_ipac",
        topic_prefix="winter",
        base_name="WNTR",
        broadcast=BROADCAST_BOOL,
        save_local=True,
        avro_schema_path=winter_avro_schema_path,
    ),
    HeaderEditor(edit_keys="sent", values=BROADCAST_BOOL),
    DatabaseSourceInserter(
        db_table=Candidate,
        duplicate_protocol="replace",
    ),
]

avro_export = avro_write + avro_broadcast

process_candidates = ml_classify + crossmatch_candidates + name_candidates + avro_write

load_avro = [SourceLoader(input_dir_name="preavro"), SourceBatcher(BASE_NAME_KEY)]

load_skyportal = [
    SourceLoader(input_dir_name="preskyportal"),
    SourceBatcher(BASE_NAME_KEY),
]

send_to_skyportal = [
    SkyportalCandidateUploader(**winter_fritz_config),
    HeaderEditor(edit_keys="sent", values=True),
    DatabaseSourceInserter(
        db_table=Candidate,
        duplicate_protocol="replace",
    ),
]

# To make a mosaic by stacking all boards
stack_boards = [
    ImageBatcher([TARGET_KEY]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=False,
        cache=False,
        center_type="ALL",
        temp_output_sub_dir="mosaic_weights",
    ),
    ImageSaver(output_dir_name="mosaic"),
]

mosaic = load_final_stack + stack_boards

# To make cals for focusing
focus_subcoord = [
    HeaderAnnotator(input_keys=["BOARD_ID"], output_key="SUBCOORD"),
    HeaderAnnotator(input_keys=["BOARD_ID"], output_key="SUBDETID"),
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

photcal_stacks = [
    ImageLoader(
        input_sub_dir="stack",
        input_img_dir=base_output_dir,
        load_image=load_winter_stack,
    ),
] + photcal_and_export

reduce_unpacked = load_and_export_unpacked + full_reduction

reduce_unpacked_subset = (
    load_unpacked + select_subset + export_unpacked + full_reduction
)

reduce = unpack_all + full_reduction

reftest = (
    unpack_subset
    + dark_calibrate
    + flat_calibrate
    + process_and_stack
    + select_ref
    + refbuild
)

detrend_unpacked = load_and_export_unpacked + dark_calibrate + flat_calibrate

only_ref = load_ref + select_ref + refbuild

realtime = extract_all + mask_and_split + save_raw + full_reduction

candidates = detect_candidates + process_candidates

full = realtime + imsub

focus_cals = (
    load_raw
    + extract_all
    + mask
    + focus_subcoord
    + csvlog
    + dark_calibrate
    + flat_calibrate
)

stack_forced_photometry = [
    ImageRebatcher([BASE_NAME_KEY]),
    ForcedPhotometryDetector(ra_header_key="TARGRA", dec_header_key="TARGDEC"),
    AperturePhotometry(
        aper_diameters=[5, 8, 10, 15],
        phot_cutout_half_size=50,
        bkg_in_diameters=[20, 20, 20, 20],
        bkg_out_diameters=[40, 40, 40, 40],
    ),
]

diff_forced_photometry = [
    ImageDebatcher(),
    ImageBatcher([BASE_NAME_KEY]),
    ForcedPhotometryDetector(ra_header_key="TARGRA", dec_header_key="TARGDEC"),
    AperturePhotometry(
        aper_diameters=[5, 8, 10, 15],
        phot_cutout_half_size=50,
        bkg_in_diameters=[20, 20, 20, 20],
        bkg_out_diameters=[40, 40, 40, 40],
    ),
    PSFPhotometry(),
]

perform_astrometry = load_calibrated + fourier_filter + astrometry
# + validate_astrometry
