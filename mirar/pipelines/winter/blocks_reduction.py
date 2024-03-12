"""
Module with blocks for WINTER image reduction and calibration
"""

# pylint: disable=duplicate-code
from mirar.paths import (
    BASE_NAME_KEY,
    DITHER_N_KEY,
    EXPTIME_KEY,
    FITS_MASK_KEY,
    LATEST_SAVE_KEY,
    MAX_DITHER_KEY,
    OBSCLASS_KEY,
    RAW_IMG_KEY,
    TARGET_KEY,
    base_output_dir,
)
from mirar.pipelines.winter.config import (
    psfex_path,
    scamp_config_path,
    sextractor_astrometry_config,
    sextractor_astromstats_config,
    sextractor_photometry_config,
    sextractor_photometry_psf_config,
    swarp_config_path,
)
from mirar.pipelines.winter.constants import NXSPLIT, NYSPLIT
from mirar.pipelines.winter.generator import (
    select_winter_flat_images,
    winter_anet_sextractor_config_path_generator,
    winter_astrometric_ref_catalog_generator,
    winter_astrometric_ref_catalog_namer,
    winter_astrometry_sextractor_catalog_purifier,
    winter_astrostat_catalog_purifier,
    winter_cal_requirements,
    winter_fourier_filtered_image_generator,
    winter_photcal_color_columns_generator,
    winter_photometric_catalog_generator,
    winter_photometric_catalogs_purifier,
    winter_photometric_ref_catalog_namer,
    winter_stackid_annotator,
)
from mirar.pipelines.winter.load_winter_image import (
    annotate_winter_subdet_headers,
    get_raw_winter_mask,
    load_winter_mef_image,
    load_winter_stack,
)
from mirar.pipelines.winter.models import (
    AstrometryStat,
    Exposure,
    FirstPassAstrometryStat,
    Raw,
    Stack,
)
from mirar.pipelines.winter.validator import (
    masked_images_rejector,
    poor_astrometric_quality_rejector,
    winter_dark_oversubtraction_rejector,
)
from mirar.processors.astromatic import PSFex, Scamp
from mirar.processors.astromatic.scamp.scamp import SCAMP_HEADER_KEY
from mirar.processors.astromatic.sextractor.background_subtractor import (
    SextractorBkgSubtractor,
)
from mirar.processors.astromatic.sextractor.sextractor import (
    Sextractor,
    sextractor_checkimg_map,
)
from mirar.processors.astromatic.swarp import ReloadSwarpComponentImages
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.astrometry.utils import AstrometryFromFile
from mirar.processors.astrometry.validate import AstrometryStatsWriter
from mirar.processors.catalog_limiting_mag import CatalogLimitingMagnitudeCalculator
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_inserter import (
    DatabaseImageBatchInserter,
    DatabaseImageInserter,
)
from mirar.processors.database.database_updater import ImageDatabaseMultiEntryUpdater
from mirar.processors.flat import FlatCalibrator, SkyFlatCalibrator
from mirar.processors.mask import (  # MaskAboveThreshold,
    MaskDatasecPixels,
    MaskPixelsFromFunction,
    MaskPixelsFromPathInverted,
    MaskPixelsFromWCS,
    WriteMaskedCoordsToFile,
)
from mirar.processors.photcal import ZPWithColorTermCalculator
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.split import SplitImage
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImagePlotter,
    ImageRejector,
    ImageSaver,
    ImageSelector,
    MEFLoader,
)
from mirar.processors.utils.cal_hunter import CalHunter
from mirar.processors.utils.image_loader import LoadImageFromHeader

# Start for new WINTER blocks:

load_raw = [
    MEFLoader(
        input_sub_dir="raw",
        load_image=load_winter_mef_image,
    ),
]

extract_all = [
    ImageBatcher("UTCTIME"),
    DatabaseImageBatchInserter(db_table=Exposure, duplicate_protocol="replace"),
    ImageSelector((OBSCLASS_KEY, ["dark", "science", "flat"])),
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
        ("BOARD_ID", str(BOARD_ID)),
    ),
]

# mask
mask = [
    ImageBatcher(BASE_NAME_KEY),
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
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
    DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace"),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", TARGET_KEY, "SUBCOORD"]),
    CustomImageBatchModifier(winter_stackid_annotator),
    ImageSaver(output_dir_name="raw_unpacked", write_mask=False),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    ImageRejector(("BOARD_ID", "0")),
]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Various processing steps

# Load from unpacked dir

load_unpacked = [
    ImageLoader(input_sub_dir="raw_unpacked", input_img_dir=base_output_dir),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    ImageRejector(("BOARD_ID", "0")),
    ImageDebatcher(),
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
    DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace"),
    ImageRejector(("BOARD_ID", "0")),
]

# Calibration hunter
cal_hunter = [
    CalHunter(load_image=load_winter_mef_image, requirements=winter_cal_requirements)
]

# Detrend blocks
dark_calibrate = [
    ImageDebatcher(),
    ImageBatcher(
        ["BOARD_ID", EXPTIME_KEY, "SUBCOORD", "GAINCOLT", "GAINCOLB", "GAINROW"]
    ),
    DarkCalibrator(
        cache_sub_dir="calibration_darks",
        cache_image_name_header_keys=[EXPTIME_KEY, "BOARD_ID"],
    ),
    ImageSelector((OBSCLASS_KEY, ["science", "flat"])),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "UTCTIME", "SUBCOORD"]),
    ImageSaver(output_dir_name="darkcal"),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    CustomImageBatchModifier(winter_dark_oversubtraction_rejector),
]

# First pass flat calibration
first_pass_flat_calibrate = [
    ImageDebatcher(),
    ImageBatcher(
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
        cache_sub_dir="fp_flats",
        select_flat_images=select_winter_flat_images,
        cache_image_name_header_keys=["FILTER", "BOARD_ID", TARGET_KEY],
    ),
    ImageSaver(output_dir_name="skyflatcal"),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "UTCTIME", "SUBCOORD"]),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="fp_skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    ImageSaver(output_dir_name="fp_skysub"),
]

# Fourier filtering
fourier_filter = [
    CustomImageBatchModifier(winter_fourier_filtered_image_generator),
]

# Various astrometry-related blocks
astrometry_net = [
    ImageDebatcher(),
    ImageBatcher(["UTCTIME", "BOARD_ID", "SUBCOORD"]),
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
        cache=True,
    ),
]

astrometry_scamp = [
    ImageDebatcher(),
    ImageBatcher(
        [TARGET_KEY, "FILTER", EXPTIME_KEY, "BOARD_ID", "SUBCOORD", "DITHGRP"]
    ),
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
        cache=False,
        crossmatch_radius_arcsec=5.0,
    ),
]

first_pass_validate_astrometry_export_and_filter = validate_astrometry + [
    DatabaseImageInserter(
        db_table=FirstPassAstrometryStat, duplicate_protocol="replace"
    ),
    CustomImageBatchModifier(poor_astrometric_quality_rejector),
]

second_pass_validate_astrometry_export_and_filter = validate_astrometry + [
    DatabaseImageInserter(db_table=AstrometryStat, duplicate_protocol="replace"),
    CustomImageBatchModifier(poor_astrometric_quality_rejector),
]

load_calibrated = [
    ImageLoader(input_sub_dir="skysub", input_img_dir=base_output_dir),
    ImageBatcher(["UTCTIME", "BOARD_ID"]),
]

# Stacking
stack_dithers = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", "EXPTIME", TARGET_KEY, "SUBCOORD"]),
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
]

first_pass_stacking = (
    astrometry_net
    + astrometry_scamp
    + first_pass_validate_astrometry_export_and_filter
    + [
        ImageSaver(output_dir_name="fp_post_astrometry"),
    ]
    + stack_dithers
    + [
        ImageSaver(output_dir_name="fp_stack"),
    ]
)

load_astrometried = [
    ImageLoader(
        input_sub_dir="fp_post_astrometry",
        input_img_dir=base_output_dir,
    )
]

# Second pass calibration
second_pass_calibration = [
    ImageLoader(
        input_sub_dir="fp_stack",
        input_img_dir=base_output_dir,
        load_image=load_winter_stack,
    ),
    ImageDebatcher(),
    ImageBatcher([TARGET_KEY, "BOARD_ID"]),
    Sextractor(
        output_sub_dir="sp_stack_source_mask",
        **sextractor_astrometry_config,
        checkimage_type="SEGMENTATION",
        cache=True,
    ),
    MaskPixelsFromPathInverted(
        mask_path_key=sextractor_checkimg_map["SEGMENTATION"],
        write_masked_pixels_to_file=True,
        output_dir="sp_stack_source_mask",
    ),
    WriteMaskedCoordsToFile(output_dir="sp_stack_mask"),
    ReloadSwarpComponentImages(
        copy_header_keys=FITS_MASK_KEY,
    ),
    LoadImageFromHeader(
        header_key=RAW_IMG_KEY,
        copy_header_keys=[SCAMP_HEADER_KEY, FITS_MASK_KEY],
    ),
    AstrometryFromFile(astrometry_file_key=SCAMP_HEADER_KEY),
    ImageSaver(output_dir_name="sp_astrometry", write_mask=True),
    MaskPixelsFromWCS(
        write_masked_pixels_to_file=True,
        output_dir="sp_source_mask",
        only_write_mask=True,
    ),
    ImageSaver(output_dir_name="sp_masked", write_mask=True),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY),
    ImageSaver(output_dir_name="sp_calibration_flat"),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    ImageSaver(output_dir_name="skysub"),
]

# Second pass astrometry
second_pass_astrometry = (
    astrometry_net
    + [
        ImageSaver(output_dir_name="post_anet"),
    ]
    + astrometry_scamp
    + [
        ImageSaver(output_dir_name="post_scamp"),
    ]
)

# Second pass stacking
second_pass_stack = (
    second_pass_astrometry
    + second_pass_validate_astrometry_export_and_filter
    + stack_dithers
    + [
        ImageSaver(output_dir_name="stack"),
    ]
)

# Photometric calibration
photcal_and_export = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD"]),
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
    ImageSaver(output_dir_name="final"),
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

photcal_stacks = [
    ImageLoader(
        input_sub_dir="stack",
        input_img_dir=base_output_dir,
        load_image=load_winter_stack,
    ),
] + photcal_and_export

# Final stack
load_final_stack = [
    ImageLoader(
        input_sub_dir="final",
        input_img_dir=base_output_dir,
        load_image=load_winter_stack,
    ),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="ignore"),
]

plot_stack = [
    ImageDebatcher(),
    ImageBatcher([TARGET_KEY, "BOARD_ID"]),
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

# End of image-reduction blocks
