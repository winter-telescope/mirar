"""
Module for WINTER data reduction
"""
from mirar.io import open_fits
from mirar.paths import (
    BASE_NAME_KEY,
    DITHER_N_KEY,
    EXPTIME_KEY,
    FITS_MASK_KEY,
    MAX_DITHER_KEY,
)
from mirar.pipelines.winter.config import (
    psfex_path,
    scamp_config_path,
    sextractor_anet_config,
    sextractor_autoastrometry_config,
    sextractor_photometry_config,
    sextractor_reference_config,
    swarp_config_path,
)
from mirar.pipelines.winter.generator import (
    winter_astrometric_catalog_generator,
    winter_astrostat_catalog_purifier,
    winter_photometric_catalog_generator,
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
from mirar.processors.astromatic import PSFex, Scamp
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.astrometry.validate import AstrometryStatsWriter
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_modifier import ModifyImageDatabaseSeqList
from mirar.processors.mask import (
    MaskAboveThreshold,
    MaskDatasecPixels,
    MaskPixelsFromFunction,
    WriteMaskedCoordsToFile,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.reference import GetReferenceImage, ProcessReference
from mirar.processors.sky import NightSkyMedianCalibrator, SkyFlatCalibrator
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
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

refbuild = [
    ImageDebatcher(),
    GetReferenceImage(
        ref_image_generator=winter_reference_generator,
    ),
    ImageSaver(output_dir_name="stacked_ref"),
]

BOARD_ID = 4
TARGET_NAME = "m39"

load_anet = [
    ImageLoader(input_sub_dir=f"anet_{BOARD_ID}", load_image=load_proc_winter_image),
    # ImageSelector(("TARGNAME", f"{TARGET_NAME}"), ("OBSTYPE", "SCIENCE")),
]

load_stack = [
    ImageLoader(input_sub_dir=f"anet_{BOARD_ID}", load_image=load_proc_winter_image),
    # ImageSelector(("TARGNAME", f"{TARGET_NAME}"), ("OBSTYPE", "SCIENCE")),
    ImageBatcher("EXPTIME"),
]

load_multiboard_stack = [
    ImageLoader(
        input_sub_dir=f"stack_all_{TARGET_NAME}", load_image=load_stacked_winter_image
    ),
    ImageSelector(
        # ("TARGNAME", f"{TARGET_NAME}"),
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

dark_cal_all_boards = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "EXPTIME", "SUBCOORD"]),
    WriteMaskedCoordsToFile(output_dir="mask_raw"),
    DarkCalibrator(cache_sub_dir="calibration"),
    ImageSaver(output_dir_name="darkcal"),
    ImageDebatcher(),
]

flat_cal = [
    # ImageSelector(("OBSTYPE", ["FOCUS", "SCIENCE", "FLAT"])),
    # ImageSelector(("TARGNAME", ["INTERESTING"])),
    # ImageBatcher(["BOARD_ID", "FILTER"]),
    # FlatCalibrator(flat_mask_key=FITS_MASK_KEY,
    #                cache_sub_dir=f"calibration_{board_id}"
    #                ),
    ImageSelector(("OBSTYPE", ["SCIENCE"]), ("TARGNAME", f"{TARGET_NAME}")),
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "EXPTIME", "SUBCOORD"]),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY, cache_sub_dir=f"skycals_{BOARD_ID}"),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageSaver(output_dir_name=f"skyflatcal_{BOARD_ID}"),
    NightSkyMedianCalibrator(flat_mask_key=FITS_MASK_KEY),
    ImageSaver(output_dir_name=f"skysub_{BOARD_ID}"),
]

flat_cal_all_boards = [
    # ImageSelector(("OBSTYPE", ["SCIENCE"]), ("TARGNAME", f"{TARGET_NAME}")),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", "SUBCOORD"]),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY, cache_sub_dir="skycals"),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageSaver(output_dir_name="skyflatcal"),
    NightSkyMedianCalibrator(flat_mask_key=FITS_MASK_KEY),
    ImageSaver(output_dir_name="skysub"),
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

process_proc_all_boards = [
    ImageDebatcher(),
    ImageBatcher(["UTCTIME", "BOARD_ID", "SUBCOORD"]),
    # ImageSelector(("FIELDID", ["2789", "0697", "9170"])),
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
    # Sextractor(
    #     **sextractor_autoastrometry_config,
    #     write_regions_bool=True,
    #     output_sub_dir="scamp",
    # ),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "SUBCOORD"]),
    # Scamp(
    #     temp_output_sub_dir="scamp",
    #     ref_catalog_generator=winter_astrometric_catalog_generator,
    #     scamp_config_path=scamp_config_path,
    #     cache=False,
    # ),
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
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "SUBCOORD"]),
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
    ImageSaver(output_dir_name="stack"),
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

stack_proc = [
    ImageBatcher(["TARGNAME", "FILTER"]),
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
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "SUBCOORD"]),
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

extract_all = [
    MEFLoader(
        input_sub_dir="raw",
        load_image=load_winter_mef_image,
        extension_num_header_key="BOARD_ID",
    ),
    ImageDebatcher(),
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
            "TARGNAME",
            "RADEG",
            "DECDEG",
            "T_ROIC",
        ]
    ),
    ImageSelector(("OBSTYPE", ["DARK", "SCIENCE"])),
    ImageBatcher(BASE_NAME_KEY),
]

select_split_subset = [ImageSelector(("SUBCOORD", "0_0"))]

make_log_and_save = []

# Optional subset selection

select_subset = [
    ImageSelector(
        ("EXPTIME", "120.0"), ("BOARD_ID", str(BOARD_ID)), ("FILTER", ["dark", "Y"])
    ),
]

# Split

mask_and_split = [
    MaskAboveThreshold(threshold=40000.0),
    MaskDatasecPixels(),
    MaskPixelsFromFunction(mask_function=get_raw_winter_mask),
    ImageDebatcher(),
    ImageBatcher("UTCTIME"),
    DatabaseImageBatchExporter(db_table=Exposures, duplicate_protocol="ignore"),
    SplitImage(n_x=NXSPLIT, n_y=NYSPLIT),
    CustomImageModifier(annotate_winter_subdet_headers),
]

# Save raw images

save_raw = [
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
    DatabaseImageExporter(db_table=Raw, duplicate_protocol="replace", q3c_bool=False),
]

unpack_subset = extract_all + select_subset + mask_and_split + save_raw
unpack_all = extract_all + mask_and_split + save_raw

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Various processing steps

# Load from unpacked dir

load_unpacked = [
    ImageLoader(input_sub_dir="raw_unpacked", load_image=open_fits),
]

# Image subtraction

imsub = [
    ImageLoader(input_sub_dir="photcal", load_image=open_fits),
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

final = [
    ImageLoader(input_sub_dir="prezogy", load_image=open_fits),
    ZOGY(output_sub_dir="subtract", sci_zp_header_key="ZP_AUTO"),
    ImageSaver(output_dir_name="diffs"),
]

full_commissioning = load_unpacked + detrend + process_detrended  # + stack_proc

full_commissioning_proc = (
    dark_cal_all_boards + flat_cal_all_boards + process_proc_all_boards  # + photcal
)

full_commissioning_all_boards = load_unpacked + full_commissioning_proc


reduce = unpack_all + full_commissioning_proc
