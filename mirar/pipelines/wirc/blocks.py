"""
Module containing standard processing blocks for WIRC
"""

# pylint: disable=duplicate-code
from mirar.paths import (
    BASE_NAME_KEY,
    FITS_MASK_KEY,
    LATEST_SAVE_KEY,
    OBSCLASS_KEY,
    RAW_IMG_KEY,
    SATURATE_KEY,
)
from mirar.pipelines.wirc.generator import (
    annotate_target_coordinates,
    wirc_astrometric_catalog_generator,
    wirc_photometric_catalog_generator,
    wirc_reference_generator,
    wirc_reference_image_resampler,
    wirc_reference_psfex,
    wirc_reference_sextractor,
    wirc_source_table_filter_annotator,
    wirc_zogy_catalogs_purifier,
)
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_files import (
    psfex_path,
    scamp_fp_path,
    sextractor_astrometry_config,
    sextractor_candidate_config,
    sextractor_photometry_config,
    sextractor_reference_config,
    swarp_sp_path,
    wirc_mask_path,
)
from mirar.processors.astromatic import Scamp, Sextractor, Swarp
from mirar.processors.astromatic.psfex import PSFex
from mirar.processors.astromatic.scamp.scamp import SCAMP_HEADER_KEY
from mirar.processors.astromatic.sextractor.sextractor import sextractor_checkimg_map
from mirar.processors.astromatic.swarp import ReloadSwarpComponentImages
from mirar.processors.astrometry.autoastrometry import AutoAstrometry
from mirar.processors.astrometry.utils import AstrometryFromFile
from mirar.processors.catalog_limiting_mag import CatalogLimitingMagnitudeCalculator
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.mask import (
    MaskAboveThreshold,
    MaskPixelsFromPath,
    MaskPixelsFromPathInverted,
    MaskPixelsFromWCS,
    WriteMaskedCoordsToFile,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import ProcessReference
from mirar.processors.sky import NightSkyMedianCalibrator
from mirar.processors.sources import (
    ForcedPhotometryDetector,
    ImageUpdater,
    JSONExporter,
    ParquetWriter,
    SourceWriter,
    ZOGYSourceDetector,
)
from mirar.processors.sources.source_table_modifier import CustomSourceTableModifier
from mirar.processors.sources.utils import RegionsWriter
from mirar.processors.utils import (
    CustomImageBatchModifier,
    HeaderAnnotator,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageRebatcher,
    ImageSaver,
    ImageSelector,
)
from mirar.processors.utils.image_loader import LoadImageFromHeader
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [ImageLoader(input_sub_dir="raw", load_image=load_raw_wirc_image)]

load_stack = [
    ImageLoader(input_sub_dir="final", load_image=load_raw_wirc_image),
    ImageBatcher(split_key=[BASE_NAME_KEY]),
]

log = [
    ImageRebatcher("UTSHUT"),
    CSVLog(
        export_keys=[
            "OBJECT",
            "FILTER",
            "UTSHUT",
            "EXPTIME",
            "COADDS",
            OBSCLASS_KEY,
            BASE_NAME_KEY,
        ]
    ),
    ImageDebatcher(),
]

masking = [
    ImageSelector((OBSCLASS_KEY, ["science", "dark"])),
    MaskPixelsFromPath(mask_path=wirc_mask_path),
]

dark_calibration = [ImageBatcher("EXPTIME"), DarkCalibrator()]


reduction = [
    ImageSaver(output_dir_name="darkcal"),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    ImageDebatcher(),
    ImageSelector((OBSCLASS_KEY, "science")),
    ImageBatcher(split_key=["filter", "object"]),
    CustomImageBatchModifier(annotate_target_coordinates),
    SkyFlatCalibrator(cache_sub_dir="firstpasscal"),
    NightSkyMedianCalibrator(cache_sub_dir="firstpasscal"),
    ImageBatcher(BASE_NAME_KEY),
    AutoAstrometry(catalog="tmc"),
    Sextractor(output_sub_dir="postprocess", **sextractor_astrometry_config),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
        cache=True,
        temp_output_sub_dir="firstpassscamp",
    ),
    ImageRebatcher(split_key=["filter", "object"]),
    ImageSaver(output_dir_name="firstpass"),
    Swarp(
        swarp_config_path=swarp_sp_path,
        calculate_dims_in_swarp=True,
        temp_output_sub_dir="firstpassswarp",
    ),
    ImageSaver(output_dir_name="firstpassstack"),
    Sextractor(
        output_sub_dir="firstpasssextractor",
        **sextractor_astrometry_config,
        checkimage_type="SEGMENTATION",
        cache=True,
    ),
    MaskPixelsFromPathInverted(
        mask_path_key=sextractor_checkimg_map["SEGMENTATION"],
        write_masked_pixels_to_file=True,
        output_dir="mask1",
    ),
    ImageSaver(output_dir_name="mask1", write_mask=True),
    MaskAboveThreshold(
        threshold_key=SATURATE_KEY, write_masked_pixels_to_file=True, output_dir="mask2"
    ),
    ImageSaver(output_dir_name="mask2", write_mask=True),
    WriteMaskedCoordsToFile(output_dir="mask_stack"),
    ReloadSwarpComponentImages(
        load_image=load_raw_wirc_image,
        copy_header_keys=[FITS_MASK_KEY, "TARGRA", "TARGDEC"],
    ),
    LoadImageFromHeader(
        header_key=RAW_IMG_KEY,
        copy_header_keys=[SCAMP_HEADER_KEY, FITS_MASK_KEY, "TARGRA", "TARGDEC"],
        load_image=load_raw_wirc_image,
    ),
    AstrometryFromFile(astrometry_file_key=SCAMP_HEADER_KEY),
    ImageSaver(output_dir_name="firstpassastrom", write_mask=True),
    MaskPixelsFromWCS(
        write_masked_pixels_to_file=True,
        output_dir="mask_secondpass",
        only_write_mask=True,
    ),
    ImageSaver(output_dir_name="firstpassmasked", write_mask=True),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY, cache_sub_dir="secondpasscal"),
    NightSkyMedianCalibrator(
        flat_mask_key=FITS_MASK_KEY, cache_sub_dir="secondpasscal"
    ),
    Sextractor(output_sub_dir="postprocess", **sextractor_astrometry_config),
    Swarp(
        swarp_config_path=swarp_sp_path,
        calculate_dims_in_swarp=True,
        temp_output_sub_dir="secondpassswarp",
    ),
    ImageSaver(output_dir_name="stack"),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir="final_sextractor",
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=wirc_photometric_catalog_generator,
        write_regions=True,
        temp_output_sub_dir="photcal",
    ),
    CatalogLimitingMagnitudeCalculator(
        sextractor_mag_key_name="MAG_AUTO", write_regions=True
    ),
    ImageSaver(output_dir_name="final"),
]

reduce = log + masking + dark_calibration + reduction

reference = [
    ImageRebatcher(split_key=[BASE_NAME_KEY]),
    ProcessReference(
        ref_image_generator=wirc_reference_generator,
        swarp_resampler=wirc_reference_image_resampler,
        sextractor=wirc_reference_sextractor,
        ref_psfex=wirc_reference_psfex,
    ),
]

subtract = [
    Sextractor(**sextractor_reference_config, output_sub_dir="subtract", cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="subtract", norm_fits=True),
    ZOGYPrepare(output_sub_dir="subtract"),
    ZOGY(output_sub_dir="subtract", catalog_purifier=wirc_zogy_catalogs_purifier),
    ImageSaver(output_dir_name="diffs"),
]

export_candidates_from_header = [
    ForcedPhotometryDetector(ra_header_key="TARGRA", dec_header_key="TARGDEC"),
    RegionsWriter(output_dir_name="diffs"),
]

calculate_photometry = [
    AperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_half_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    PSFPhotometry(),
]

export_photometry = [
    SourceWriter(output_dir_name="photometry"),
    JSONExporter(output_dir_name="photometry"),
    ParquetWriter(output_dir_name="photometry"),
    ImageUpdater(modify_dir_name="diffs"),
]

forced_photometry = (
    export_candidates_from_header + calculate_photometry + export_photometry
)

image_photometry = [
    Sextractor(**sextractor_reference_config, output_sub_dir="subtract", cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="photometry", norm_fits=True),
    ImageSaver(output_dir_name="photometry"),
]

imsub = reference + subtract + forced_photometry

candidates = (
    [
        ImageLoader(input_sub_dir="diffs"),
        ZOGYSourceDetector(
            output_sub_dir="subtract",
            **sextractor_candidate_config,
        ),
        RegionsWriter(output_dir_name="candidates"),
    ]
    + calculate_photometry
    + [
        CustomSourceTableModifier(modifier_function=wirc_source_table_filter_annotator),
        JSONExporter(output_dir_name="candidates"),
        ParquetWriter(output_dir_name="candidates"),
        SourceWriter(output_dir_name="candidates"),
    ]
)

test = reduction + imsub
