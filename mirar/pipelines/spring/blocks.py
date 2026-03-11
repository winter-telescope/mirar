from mirar.paths import (
    BASE_NAME_KEY,
    EXPTIME_KEY,
    LATEST_SAVE_KEY,
    OBSCLASS_KEY,
    RAW_IMG_KEY,
    TARGET_KEY,
    ZP_KEY,
    base_output_dir,
)
from mirar.pipelines.spring.config import (
    psfex_config_path,
    ref_psfex_path,
    sextractor_astrometry_config,
    sextractor_photometry_config,
    sextractor_PSF_photometry_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
)
from mirar.pipelines.spring.generator import (
    mask_stamps_around_bright_stars,
    spring_anet_sextractor_config_path_generator,
    spring_imsub_catalog_purifier,
    spring_photcal_color_columns_generator,
    spring_photometric_catalog_generator,
    spring_ref_photometric_catalogs_purifier,
    spring_reference_generator,
    spring_reference_image_resampler_for_zogy,
    spring_reference_psf_phot_sextractor,
    spring_reference_psfex,
    spring_reference_sextractor,
    spring_stackid_annotator,
)
from mirar.pipelines.spring.load_spring_image import (
    load_raw_spring_image,
    load_spring_stack,
)
from mirar.pipelines.spring.models import Diff, Raw, Stack
from mirar.processors.astromatic import PSFex
from mirar.processors.astromatic.sextractor.background_subtractor import (
    SextractorBkgSubtractor,
)
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.catalog_limiting_mag import CatalogLimitingMagnitudeCalculator
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database import DatabaseImageInserter
from mirar.processors.database.database_updater import ImageDatabaseMultiEntryUpdater
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.mask import MaskPixelsFromFunction
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.photcal.zp_calculator import (
    OutlierRejectionZPCalculator,
    ZPWithColorTermCalculator,
)
from mirar.processors.reference import ProcessReference
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    ImageLoader,
    ImageRebatcher,
    ImageSaver,
    ImageSelector,
    ModeMasker,
)
from mirar.processors.zogy.reference_aligner import AlignReference
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [
    ImageLoader(input_sub_dir="raw", load_image=load_raw_spring_image),
    ImageRebatcher(
        [
            "FILTER",
            "EXPTIME",
            TARGET_KEY,
        ]
    ),
    CustomImageBatchModifier(spring_stackid_annotator),
    ImageSaver(output_dir_name="loaded_raw"),
    ImageRebatcher(BASE_NAME_KEY),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    DatabaseImageInserter(db_table=Raw, duplicate_protocol="replace"),
]

csvlog = [
    ImageRebatcher(BASE_NAME_KEY),
    CSVLog(
        export_keys=[
            # Time / Identity
            "UTCISO",
            "DATE-OBS",
            "MJD-OBS",
            BASE_NAME_KEY,
            "FILENAME",
            # Targeting
            TARGET_KEY,
            "TARGNAME",
            "FIELDID",
            # Pointing
            "RADEG",
            "DECDEG",
            "AZIMUTH",
            "ALTITUDE",
            "AIRMASS",
            # Instrument / exposure
            "INSTRUME",
            "TELESCOP",
            "OBSERVAT",
            "FILTER",
            EXPTIME_KEY,
            OBSCLASS_KEY,
            # Other keys
            "FOCPOS",
            "TMP_CUR",
            "TMP_SET",
            "READOUTM",
            "MEDCOUNT",
        ]
    ),
]

dark_calibrate = [
    ImageRebatcher([EXPTIME_KEY]),
    DarkCalibrator(
        cache_sub_dir="calibration_darks",
        cache_image_name_header_keys=[EXPTIME_KEY],
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="darkcal"),
    ImageSelector((OBSCLASS_KEY, ["science"])),
]

flat_calibrate = [
    ImageRebatcher([TARGET_KEY]),
    SkyFlatCalibrator(
        cache_sub_dir="sky_dither_flats",
        cache_image_name_header_keys=["FILTER", TARGET_KEY],
        flat_mode="median",
        try_load_cache=False,
    ),
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

astrometry = [
    ImageRebatcher(BASE_NAME_KEY),
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[0.3, 0.5],
        scale_units="app",
        use_sextractor=True,
        # parity="neg",
        search_radius_deg=0.15,
        sextractor_config_path=spring_anet_sextractor_config_path_generator,
        use_weight=True,
        timeout=120,
        cache=False,
        no_tweak=True,
    ),
    ImageSaver("post_astrometry"),
]

stack_dithers = [
    ImageRebatcher([TARGET_KEY]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=False,
        center_type="MOST",
        temp_output_sub_dir="stacks_weights",
        header_keys_to_combine=["RAWID"],
        min_required_coadds=3,
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ModeMasker(),
    ImageSaver(output_dir_name="stack"),
]

photcal_with_color = [
    ImageRebatcher(BASE_NAME_KEY),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="stack_psf",
        checkimage_name="BACKGROUND_RMS",
    ),
    PSFex(
        config_path=psfex_config_path,
        output_sub_dir="phot",
        norm_fits=True,
    ),
    Sextractor(
        **sextractor_PSF_photometry_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
        use_psfex=True,
    ),
    PhotCalibrator(
        ref_catalog_generator=spring_photometric_catalog_generator,
        catalogs_purifier=spring_ref_photometric_catalogs_purifier,
        zp_calculator=ZPWithColorTermCalculator(
            color_colnames_guess_generator=spring_photcal_color_columns_generator,
            reject_outliers=True,
            solver="curve_fit",
        ),
        write_regions=True,
        zp_column_name="MAG_POINTSOURCE",
    ),
    CatalogLimitingMagnitudeCalculator(
        sextractor_mag_key_name="MAG_POINTSOURCE", write_regions=True
    ),
    ImageSaver(output_dir_name="processed_after_psf_with_color"),
]

photcal_without_color = [
    ImageRebatcher(BASE_NAME_KEY),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="stack_psf",
        checkimage_name="BACKGROUND_RMS",
    ),
    PSFex(
        config_path=psfex_config_path,
        output_sub_dir="phot",
        norm_fits=True,
    ),
    Sextractor(
        **sextractor_PSF_photometry_config,
        output_sub_dir="phot",
        checkimage_type="BACKGROUND_RMS",
        use_psfex=True,
    ),
    PhotCalibrator(
        ref_catalog_generator=spring_photometric_catalog_generator,
        catalogs_purifier=spring_ref_photometric_catalogs_purifier,
        zp_calculator=OutlierRejectionZPCalculator(),
        write_regions=True,
        zp_column_name="MAG_POINTSOURCE",
    ),
    CatalogLimitingMagnitudeCalculator(
        sextractor_mag_key_name="MAG_AUTO", write_regions=True
    ),
    ImageSaver(output_dir_name="processed_after_psf"),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="replace"),
    ImageDatabaseMultiEntryUpdater(
        sequence_key="rawid",
        db_table=Raw,
        db_alter_columns="ustackid",
    ),
]

load_final_stack = [
    ImageLoader(
        input_sub_dir="processed_after_psf",
        input_img_dir=base_output_dir,
        load_image=load_spring_stack,
    ),
    DatabaseImageInserter(db_table=Stack, duplicate_protocol="ignore"),
    ImageRebatcher(BASE_NAME_KEY),
]

imsub = [
    ImageRebatcher(["FILTER", TARGET_KEY, "STACKID"]),
    ProcessReference(
        ref_image_generator=spring_reference_generator,
        swarp_resampler=spring_reference_image_resampler_for_zogy,
        sextractor=spring_reference_sextractor,
        ref_psfex=spring_reference_psfex,
        phot_sextractor=spring_reference_psf_phot_sextractor,
    ),
    Sextractor(
        **sextractor_reference_psf_phot_config,
        output_sub_dir="subtract",
        cache=False,
        use_psfex=True,
    ),
    PSFex(config_path=ref_psfex_path, output_sub_dir="subtract", norm_fits=True),
    AlignReference(
        order=1,
        sextractor=spring_reference_sextractor,
        psfex=spring_reference_psfex,
        phot_sextractor=spring_reference_psf_phot_sextractor,
        catalog_purifier=spring_imsub_catalog_purifier,
    ),
    MaskPixelsFromFunction(mask_function=mask_stamps_around_bright_stars),
    ImageSaver(output_dir_name="presubtract"),
    ZOGYPrepare(
        output_sub_dir="subtract",
        sci_zp_header_key="ZP_AUTO",
        ref_zp_header_key=ZP_KEY,
        catalog_purifier=spring_imsub_catalog_purifier,
        x_key="XMODEL_IMAGE",
        y_key="YMODEL_IMAGE",
        flux_key="FLUX_POINTSOURCE",
    ),
    ZOGY(
        output_sub_dir="subtract",
        sci_zp_header_key="ZP_AUTO",
        ref_zp_header_key=ZP_KEY,
    ),
    ImageSaver(output_dir_name="diffs"),
    DatabaseImageInserter(db_table=Diff, duplicate_protocol="replace"),
    ImageSaver(output_dir_name="subtract"),
]
