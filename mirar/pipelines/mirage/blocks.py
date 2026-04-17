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
from mirar.pipelines.mirage.config import (
    MIRAGE_GAIN,
    psfex_config_path,
    ref_psfex_path,
    scamp_path,
    sextractor_astrometry_config,
    sextractor_candidates_config,
    sextractor_photometry_config,
    sextractor_PSF_photometry_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
)
from mirar.pipelines.mirage.generator import (
    mask_stamps_around_bright_stars,
    mirage_anet_sextractor_config_path_generator,
    mirage_astrometric_ref_catalog_generator,
    mirage_astrometry_sextractor_catalog_purifier,
    mirage_candidate_annotator_filterer,
    mirage_imsub_catalog_purifier,
    mirage_masterbias_path_generator,
    mirage_masterdark_path_generator,
    mirage_photcal_color_columns_generator,
    mirage_photometric_catalog_generator,
    mirage_ref_photometric_catalogs_purifier,
    mirage_reference_generator,
    mirage_reference_image_resampler_for_zogy,
    mirage_reference_psf_phot_sextractor,
    mirage_reference_psfex,
    mirage_reference_sextractor,
    mirage_stack_gain_modifier,
    mirage_stackid_annotator,
)
from mirar.pipelines.mirage.load_mirage_image import (
    load_mirage_cal_image,
    load_mirage_stack,
    load_raw_mirage_image,
)
from mirar.processors.astromatic import PSFex, Scamp
from mirar.processors.astromatic.sextractor.background_subtractor import (
    SextractorBkgSubtractor,
)
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet import AstrometryNet
from mirar.processors.bias import MasterBiasCalibrator
from mirar.processors.catalog_limiting_mag import CatalogLimitingMagnitudeCalculator
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import MasterDarkCalibrator
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.mask import MaskPixelsFromFunction
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.photcal.zp_calculator import (
    OutlierRejectionZPCalculator,
    ZPWithColorTermCalculator,
)
from mirar.processors.photometry.aperture_photometry import AperturePhotometry
from mirar.processors.photometry.psf_photometry import PSFPhotometry
from mirar.processors.reference import ProcessReference
from mirar.processors.sources.forced_photometry import ForcedPhotometryDetector
from mirar.processors.sources.source_detector import ZOGYSourceDetector
from mirar.processors.sources.source_exporter import SourceWriter
from mirar.processors.sources.source_loader import SourceLoader
from mirar.processors.sources.source_selector import SourceBatcher
from mirar.processors.sources.source_table_modifier import CustomSourceTableModifier
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageRebatcher,
    ImageSaver,
    ImageSelector,
    ModeMasker,
)
from mirar.processors.utils.image_plotter import ImagePlotter
from mirar.processors.zogy.reference_aligner import AlignReference
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [
    ImageLoader(input_sub_dir="raw", load_image=load_raw_mirage_image),
    ImageRebatcher(BASE_NAME_KEY),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
]


process_cals = [
    ImageLoader(input_sub_dir="cals", load_image=load_mirage_cal_image),
    ImageRebatcher(BASE_NAME_KEY),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    ImageSaver("mirared_cals"),
]

csvlog = [
    ImageRebatcher(BASE_NAME_KEY),
    CSVLog(
        export_keys=[
            # Time / Identity
            "DATE-OBS",
            "MJD-OBS",
            BASE_NAME_KEY,
            # Targeting
            TARGET_KEY,
            "RADEG",
            "DECDEG",
            # Instrument / exposure
            "INSTRUME",
            "FILTER",
            EXPTIME_KEY,
            OBSCLASS_KEY,
            # Other keys
            "MEDCOUNT",
        ]
    ),
]

dark_calibrate = [
    ImageRebatcher([EXPTIME_KEY]),
    MasterBiasCalibrator(master_image_path_generator=mirage_masterbias_path_generator),
    MasterDarkCalibrator(
        master_image_path_generator=mirage_masterdark_path_generator,
    ),
    ImageRebatcher(BASE_NAME_KEY),
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
        scale_bounds=[0.15, 0.4],
        scale_units="app",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=2,
        sextractor_config_path=mirage_anet_sextractor_config_path_generator,
        use_weight=True,
        timeout=120,
        cache=False,
        no_tweak=True,
    ),
]

# astrometry = [
#     ImageRebatcher(BASE_NAME_KEY),
#     Sextractor(
#         **sextractor_astrometry_config,
#         write_regions_bool=True,
#         output_sub_dir="scamp",
#         cache=False,
#         catalog_purifier=mirage_astrometry_sextractor_catalog_purifier,
#         verbose_type="FULL",
#     ),
#     ImageRebatcher([TARGET_KEY, "FILTER", EXPTIME_KEY]),
#     Scamp(
#         scamp_config_path=scamp_path,
#         ref_catalog_generator=mirage_astrometric_ref_catalog_generator,
#         copy_scamp_header_to_image=True,
#         cache=False,
#         make_checkplots=True,
#     ),
#     # AutoAstrometry(catalog="tmc"),
#     # ImageSaver("post_astrometry"),
# ]

stack_dithers = [
    ImageRebatcher([TARGET_KEY, "FILTER"]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=False,
        center_type="MOST",
        temp_output_sub_dir="stacks_weights",
        # header_keys_to_combine=["RAWID"],
        min_required_coadds=3,
        gain=MIRAGE_GAIN,
    ),
    CustomImageBatchModifier(mirage_stack_gain_modifier),
    ImageRebatcher(BASE_NAME_KEY),
    ModeMasker(),
    ImageSaver(output_dir_name="stack"),
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
        ref_catalog_generator=mirage_photometric_catalog_generator,
        catalogs_purifier=mirage_ref_photometric_catalogs_purifier,
        zp_calculator=OutlierRejectionZPCalculator(),
        write_regions=True,
        zp_column_name="MAG_AUTO",
        num_matches_threshold=3,
        crossmatch_radius_arcsec=2.0,
    ),
    CatalogLimitingMagnitudeCalculator(
        sextractor_mag_key_name="MAG_AUTO", write_regions=True
    ),
    ImageSaver(output_dir_name="photcaled_stacks"),
]

export_stacks = [
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

photcal_and_export = photcal_without_color + export_stacks

load_stack = [
    ImageLoader(
        input_sub_dir="stack",
        input_img_dir=base_output_dir,
        load_image=load_mirage_stack,
    ),
    ImageRebatcher(BASE_NAME_KEY),
]

load_final_stack = [
    ImageLoader(
        input_sub_dir="processed_after_psf",
        input_img_dir=base_output_dir,
        load_image=load_mirage_stack,
    ),
    ImageRebatcher(BASE_NAME_KEY),
]

imsub = [
    ImageRebatcher(["FILTER", TARGET_KEY, "STACKID"]),
    ProcessReference(
        ref_image_generator=mirage_reference_generator,
        swarp_resampler=mirage_reference_image_resampler_for_zogy,
        sextractor=mirage_reference_sextractor,
        ref_psfex=mirage_reference_psfex,
        phot_sextractor=mirage_reference_psf_phot_sextractor,
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
        sextractor=mirage_reference_sextractor,
        psfex=mirage_reference_psfex,
        phot_sextractor=mirage_reference_psf_phot_sextractor,
        catalog_purifier=mirage_imsub_catalog_purifier,
    ),
    MaskPixelsFromFunction(mask_function=mask_stamps_around_bright_stars),
    ImageSaver(output_dir_name="presubtract"),
    ZOGYPrepare(
        output_sub_dir="subtract",
        sci_zp_header_key="ZP_AUTO",
        ref_zp_header_key=ZP_KEY,
        catalog_purifier=mirage_imsub_catalog_purifier,
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
    ImageSaver(output_dir_name="subtract"),
]

load_sub = [
    ImageLoader(input_sub_dir="diffs", input_img_dir=base_output_dir),
    ImageBatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="subtract"),
]

detect_candidates = [
    ZOGYSourceDetector(
        output_sub_dir="subtract",
        **sextractor_candidates_config,
        write_regions=True,
        detect_negative_sources=False,
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
    CustomSourceTableModifier(mirage_candidate_annotator_filterer),
    SourceWriter(output_dir_name="candidates"),
]

load_sources = [
    SourceLoader(input_dir_name="candidates", input_dir=base_output_dir),
    SourceBatcher(BASE_NAME_KEY),
]

stack_forced_photometry = [
    ImageRebatcher([BASE_NAME_KEY]),
    ForcedPhotometryDetector(ra_header_key="TARGRA", dec_header_key="TARGDEC"),
    AperturePhotometry(
        aper_diameters=[5, 8, 10, 15],
        phot_cutout_half_size=50,
        bkg_in_diameters=[20, 20, 20, 20],
        bkg_out_diameters=[40, 40, 40, 40],
    ),
    SourceWriter(output_dir_name="photometry"),
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
    SourceWriter(output_dir_name="photometry"),
]

reduce = load_raw + csvlog + dark_calibrate + flat_calibrate
