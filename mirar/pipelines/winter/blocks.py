"""
Module for WINTER data reduction
"""
from mirar.catalog.kowalski import PS1, TMASS
from mirar.downloader.get_test_data import get_test_data_dir
from mirar.io import open_raw_image
from mirar.paths import (
    BASE_NAME_KEY,
    DITHER_N_KEY,
    EXPTIME_KEY,
    FITS_MASK_KEY,
    MAX_DITHER_KEY,
    TARGET_KEY,
    base_output_dir,
    get_output_dir,
)
from mirar.pipelines.winter.config import (
    psfex_path,
    scamp_config_path,
    sextractor_anet_config,
    sextractor_autoastrometry_config,
    sextractor_candidate_config,
    sextractor_photometry_config,
    sextractor_reference_config,
    swarp_config_path,
    winter_candidate_config,
)
from mirar.pipelines.winter.generator import (
    winter_astrometric_catalog_generator,
    winter_astrostat_catalog_purifier,
    winter_photometric_catalog_generator,
    winter_refbuild_reference_generator,
    winter_reference_generator,
    winter_reference_image_resampler,
    winter_reference_psfex,
    winter_reference_sextractor,
    winter_stackid_annotator,
)
from mirar.pipelines.winter.load_winter_image import (
    annotate_winter_subdet_headers,
    get_raw_winter_mask,
    load_proc_winter_image,
    load_stacked_winter_image,
    load_winter_mef_image,
)
from mirar.pipelines.winter.models import (
    NXSPLIT,
    NYSPLIT,
    AstrometryStats,
    Exposures,
    Raw,
    Stacks,
)
from mirar.processors.alerts import AvroPacketMaker
from mirar.processors.astromatic import PSFex, Scamp
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.astrometry.validate import AstrometryStatsWriter
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_exporter import DatabaseDataframeExporter
from mirar.processors.database.database_importer import DatabaseHistoryImporter
from mirar.processors.database.database_modifier import ModifyImageDatabaseSeqList
from mirar.processors.database.utils import get_column_names_from_schema
from mirar.processors.mask import (
    MaskAboveThreshold,
    MaskDatasecPixels,
    MaskPixelsFromFunction,
    WriteMaskedCoordsToFile,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.photometry.aperture_photometry import CandidateAperturePhotometry
from mirar.processors.photometry.psf_photometry import CandidatePSFPhotometry
from mirar.processors.reference import GetReferenceImage, ProcessReference
from mirar.processors.sky import NightSkyMedianCalibrator, SkyFlatCalibrator
from mirar.processors.sources import CandidateNamer, DataframeWriter, SourceDetector
from mirar.processors.split import SUB_ID_KEY, SplitImage
from mirar.processors.sqldatabase.database_exporter import (
    DatabaseImageBatchExporter,
    DatabaseImageExporter,
)
from mirar.processors.utils import (
    CustomImageModifier,
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
    ),
    ImageSaver("testdata", output_dir=get_test_data_dir()),
]

load_test = [
    ImageLoader(
        input_img_dir=get_test_data_dir(),
        input_sub_dir="raw",
        load_image=open_raw_image,
    ),
]

refbuild = [
    GetReferenceImage(ref_image_generator=winter_refbuild_reference_generator),
    ImageSaver(output_dir_name="stacked_ref"),
]

BOARD_ID = 4
TARGET_NAME = "m39"

# Everything in this block is deprecated/untested/probably outdated

load_anet = [
    ImageLoader(input_sub_dir=f"anet_{BOARD_ID}", load_image=load_proc_winter_image),
    # ImageSelector((TARGET_KEY, f"{TARGET_NAME}"), ("OBSTYPE", "SCIENCE")),
]

load_ref = [
    ImageLoader(
        input_sub_dir="stack",
        load_image=load_stacked_winter_image,
        input_img_dir=base_output_dir,
    )
]

select_ref = [
    ImageSelector(
        ("FIELDID", str(3944)),
        ("BOARD_ID", str(BOARD_ID)),
    ),
    ImageDebatcher(),
    ImageBatcher("STACKID"),
]

load_multiboard_stack = [
    ImageLoader(
        input_sub_dir=f"stack_all_{TARGET_NAME}", load_image=load_stacked_winter_image
    ),
    ImageSelector(
        # (TARGET_KEY, f"{TARGET_NAME}"),
        ("OBSTYPE", "SCIENCE"),
    ),
]

dark_cal = [
    ImageSelector(("BOARD_ID", f"{BOARD_ID}")),
    ImageBatcher(["BOARD_ID", "EXPTIME", "SUBCOORD"]),
    WriteMaskedCoordsToFile(output_dir="mask_raw"),
    DarkCalibrator(cache_sub_dir=f"calibration_{BOARD_ID}"),
    ImageSaver(output_dir_name=f"darkcal_{BOARD_ID}"),
    ImageDebatcher(),
]

flat_cal = [
    # ImageSelector(("OBSTYPE", ["FOCUS", "SCIENCE", "FLAT"])),
    # ImageSelector((TARGET_KEY, ["INTERESTING"])),
    # ImageBatcher(["BOARD_ID", "FILTER"]),
    # FlatCalibrator(flat_mask_key=FITS_MASK_KEY,
    #                cache_sub_dir=f"calibration_{board_id}"
    #                ),
    ImageSelector(("OBSTYPE", ["SCIENCE"]), (TARGET_KEY, f"{TARGET_NAME}")),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "EXPTIME", "SUBCOORD"]),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY, cache_sub_dir=f"skycals_{BOARD_ID}"),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageSaver(output_dir_name=f"skyflatcal_{BOARD_ID}"),
    NightSkyMedianCalibrator(flat_mask_key=FITS_MASK_KEY),
    ImageSaver(output_dir_name=f"skysub_{BOARD_ID}"),
]

detrend = dark_cal + flat_cal
process_detrended = [
    ImageDebatcher(),
    AstrometryNet(
        output_sub_dir=f"anet_{BOARD_ID}",
        scale_bounds=[25, 40],
        scale_units="amw",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=1.0,
        sextractor_config_path=sextractor_autoastrometry_config["config_path"],
        use_weight=True,
        timeout=30,
    ),
    ImageSaver(output_dir_name=f"anet_{BOARD_ID}"),
    Sextractor(
        **sextractor_autoastrometry_config,
        write_regions_bool=True,
        output_sub_dir="scamp",
    ),
    Scamp(
        temp_output_sub_dir="scamp",
        ref_catalog_generator=winter_astrometric_catalog_generator,
        scamp_config_path=scamp_config_path,
        cache=False,
    ),
    ImageDebatcher(),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=True,
        subtract_bkg=False,
        cache=True,
        center_type="ALL",
        temp_output_sub_dir=f"stack_all_{TARGET_NAME}",
    ),
]

export_proc = [
    DatabaseImageExporter(db_table=Stacks, duplicate_protocol="replace", q3c_bool=False)
]

stack_proc = [
    ImageBatcher([TARGET_KEY, "FILTER"]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=True,
        cache=True,
        center_type="MANUAL",
    ),
    ImageSaver(output_dir_name=f"stack_{TARGET_NAME}"),
]

photcal = [
    # ImageSelector(("BOARD_ID", board_id)),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        temp_output_sub_dir="phot",
        write_regions=True,
        cache=True,
    ),
    AstrometryStatsWriter(
        ref_catalog_generator=winter_photometric_catalog_generator,
        write_regions=True,
        cache=True,
    ),
    # DatabaseImageExporter(
    #     db_table=AstrometryStats, duplicate_protocol="replace", q3c_bool=False
    # ),
    # ImageSaver(output_dir_name=f"phot_{board_id}_{target_name}")
    ImageSaver(output_dir_name="photcal"),
]

photcal_indiv = [
    ImageSelector(("BOARD_ID", BOARD_ID)),
    ImageDebatcher(),
    ImageBatcher(["UTCTIME"]),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir=f"phot_{BOARD_ID}_{TARGET_NAME}",
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        temp_output_sub_dir=f"phot_{BOARD_ID}_{TARGET_NAME}",
        write_regions=True,
        cache=True,
    ),
    ImageSaver(output_dir_name=f"phot_{BOARD_ID}_{TARGET_NAME}"),
]

stack_multiboard = [
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=False,
        cache=True,
        temp_output_sub_dir=f"multiboard_stack_{TARGET_NAME}",
        center_type="MANUAL",
    )
]

commissioning_multiboard_stack = load_multiboard_stack + stack_multiboard
commissioning_photcal = load_multiboard_stack + photcal
commissioning_photcal_indiv = load_anet + photcal_indiv

# Start for new WINTER blocks:

# Loading

load_raw = [
    MEFLoader(
        input_sub_dir="raw",
        load_image=load_winter_mef_image,
    ),
]

extract_all = [
    ImageSelector(("OBSTYPE", ["DARK", "SCIENCE"])),
    ImageBatcher("UTCTIME"),
    DatabaseImageBatchExporter(db_table=Exposures, duplicate_protocol="ignore"),
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
            "OBSTYPE",
            "BOARD_ID",
            "OBSCLASS",
            "TARGET",
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

select_subset = [
    ImageSelector(
        ("EXPTIME", "120.0"),
        ("FIELDID", ["3944", "999999999"]),
        ("BOARD_ID", str(BOARD_ID)),
        ("FILTER", ["dark", "J"]),
    ),
]

# Split

mask_and_split = [
    ImageBatcher(BASE_NAME_KEY),
    MaskAboveThreshold(threshold=40000.0),
    MaskDatasecPixels(),
    MaskPixelsFromFunction(mask_function=get_raw_winter_mask),
    SplitImage(n_x=NXSPLIT, n_y=NYSPLIT),
    CustomImageModifier(annotate_winter_subdet_headers),
]

# Save raw images

save_raw = [
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
    DatabaseImageExporter(db_table=Raw, duplicate_protocol="replace", q3c_bool=False),
]

unpack_subset = (
    load_raw + extract_all + csvlog + select_subset + mask_and_split + save_raw
)

unpack_all = load_raw + extract_all + csvlog + mask_and_split + save_raw

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Various processing steps

# Load from unpacked dir

load_unpacked = [
    ImageLoader(input_sub_dir="raw_unpacked"),
]

# Detrend blocks

dark_cal_all_boards = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "EXPTIME", "SUBCOORD"]),
    DarkCalibrator(cache_sub_dir="calibration"),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageSaver(output_dir_name="darkcal"),
    ImageDebatcher(),
]

flat_cal_all_boards = [
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", "SUBCOORD"]),
    SkyFlatCalibrator(cache_sub_dir="skycals"),
    ImageSaver(output_dir_name="skyflatcal"),
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", "SUBCOORD", TARGET_KEY]),
    NightSkyMedianCalibrator(),
    ImageSaver(output_dir_name="skysub"),
]

process_stack_all_boards = [
    ImageDebatcher(),
    ImageBatcher(["UTCTIME", "BOARD_ID", "SUBCOORD", "EXPTIME"]),
    ImageSaver(output_dir_name="pre_anet"),
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[15, 23],
        scale_units="amw",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=1.0,
        sextractor_config_path=sextractor_anet_config["config_path"],
        use_weight=True,
    ),
    ImageSaver(output_dir_name="post_anet"),
    Sextractor(
        **sextractor_autoastrometry_config,
        write_regions_bool=True,
        output_sub_dir="scamp",
    ),
    AstrometryStatsWriter(
        ref_catalog_generator=winter_photometric_catalog_generator,
        image_catalog_purifier=winter_astrostat_catalog_purifier,
        write_regions=True,
        cache=True,
        crossmatch_radius_arcsec=5.0,
    ),
    DatabaseImageExporter(db_table=AstrometryStats, duplicate_protocol="ignore"),
    ImageSaver(output_dir_name="anet"),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
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
    CustomImageModifier(winter_stackid_annotator),
    ImageSaver(output_dir_name="stack"),
]

photcal_and_export = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        temp_output_sub_dir="phot",
        write_regions=True,
        cache=True,
    ),
    ImageSaver(output_dir_name="final"),
    DatabaseImageExporter(
        db_table=Stacks, duplicate_protocol="replace", q3c_bool=False
    ),
    ModifyImageDatabaseSeqList(
        db_name="winter",
        schema_path="fake_placeholder_path.sql",
        sequence_key="rawid",
        db_table=Raw.sql_model.__tablename__,
        db_alter_columns="ustackid",
    ),
]

# Image subtraction

load_stack = [
    ImageLoader(input_sub_dir="final"),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
]

imsub = [
    HeaderAnnotator(input_keys=[SUB_ID_KEY], output_key="SUBDETID"),
    ProcessReference(
        ref_image_generator=winter_reference_generator,
        swarp_resampler=winter_reference_image_resampler,
        sextractor=winter_reference_sextractor,
        ref_psfex=winter_reference_psfex,
    ),
    Sextractor(**sextractor_reference_config, output_sub_dir="subtract", cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="subtract", norm_fits=True),
    ImageSaver(output_dir_name="presubtract"),
    ZOGYPrepare(output_sub_dir="subtract", sci_zp_header_key="ZP_AUTO"),
    ImageSaver(output_dir_name="prezogy"),
    ZOGY(output_sub_dir="subtract", sci_zp_header_key="ZP_AUTO"),
    ImageSaver(output_dir_name="diffs"),
]

detect_candidates = [
    HeaderAnnotator(input_keys=["ZP_AUTO"], output_key="ZP"),
    HeaderAnnotator(input_keys=["ZP_AUTO_STD"], output_key="ZP_STD"),
    SourceDetector(output_sub_dir="subtract", **sextractor_candidate_config),
]

candidate_colnames = get_column_names_from_schema(winter_candidate_config)

process_candidates = [
    DataframeWriter(output_dir_name="candidates"),
    CandidatePSFPhotometry(
        zp_colname="ZP",
    ),
    CandidateAperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
        zp_colname="ZP",
    ),
    DataframeWriter(output_dir_name="candidates"),
    XMatch(catalog=TMASS(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1(num_sources=3, search_radius_arcmin=0.5)),
    DataframeWriter(output_dir_name="kowalski"),
    DatabaseHistoryImporter(
        crossmatch_radius_arcsec=2.0,
        time_field_name="jd",
        history_duration_days=500.0,
        db_name="winter",
        db_table="candidates",
        db_output_columns=candidate_colnames,
        schema_path=winter_candidate_config,
        q3c_bool=False,
    ),
    CandidateNamer(
        db_name="winter",
        db_table="candidates",
        base_name="WNTR",
        name_start="aaaaa",
        xmatch_radius_arcsec=2,
        schema_path=winter_candidate_config,
    ),
    DatabaseDataframeExporter(
        db_name="winter",
        db_table="candidates",
        schema_path=winter_candidate_config,
        duplicate_protocol="replace",
    ),
    # DataframeWriter(output_dir_name="dbop"),
]

package_candidates = [
    AvroPacketMaker(
        output_sub_dir="avro", base_name="WNTR", broadcast=False, save_local=True
    ),
]

candidates = detect_candidates + process_candidates + package_candidates

full_commissioning = load_unpacked + detrend + process_detrended  # + stack_proc

full_commissioning_proc = (
    dark_cal_all_boards
    + flat_cal_all_boards
    + process_stack_all_boards
    + photcal_and_export
)

full_commissioning_all_boards = load_unpacked + full_commissioning_proc

reduce = unpack_all + full_commissioning_proc

reftest = (
    unpack_subset
    + dark_cal_all_boards
    + flat_cal_all_boards
    + process_stack_all_boards
    + select_ref
    + refbuild
)

only_ref = load_ref + select_ref + refbuild

realtime = extract_all + mask_and_split + save_raw + full_commissioning_proc

full = realtime + imsub
