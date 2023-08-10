"""
Module containing standard processing blocks for WIRC
"""
# pylint: disable=duplicate-code
from mirar.catalog.kowalski import PS1, TMASS
from mirar.paths import (
    CAND_NAME_KEY,
    FITS_MASK_KEY,
    LATEST_SAVE_KEY,
    OBSCLASS_KEY,
    RAW_IMG_KEY,
    SATURATE_KEY,
)
from mirar.pipelines.wirc.generator import (
    wirc_astrometric_catalog_generator,
    wirc_photometric_catalog_generator,
    wirc_photometric_img_catalog_purifier,
    wirc_reference_image_generator,
    wirc_reference_image_resampler,
    wirc_reference_psfex,
    wirc_reference_sextractor,
    wirc_source_table_filter_annotator,
    wirc_zogy_catalogs_purifier,
)
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_files import (
    prv_candidate_cols,
    psfex_path,
    scamp_fp_path,
    sextractor_astrometry_config,
    sextractor_candidate_config,
    sextractor_photometry_config,
    sextractor_reference_config,
    swarp_sp_path,
    wirc_avro_schema_path,
    wirc_mask_path,
)
from mirar.pipelines.wirc.wirc_files.models import (
    CANDIDATE_PREFIX,
    NAME_START,
    Candidate,
)
from mirar.processors.astromatic import Scamp, Sextractor, Swarp
from mirar.processors.astromatic.psfex import PSFex
from mirar.processors.astromatic.scamp.scamp import SCAMP_HEADER_KEY
from mirar.processors.astromatic.sextractor.sextractor import sextractor_checkimg_map
from mirar.processors.astromatic.swarp import ReloadSwarpComponentImages
from mirar.processors.astrometry.autoastrometry import AutoAstrometry
from mirar.processors.astrometry.utils import AstrometryFromFile
from mirar.processors.avro import IPACAvroExporter, SendToFritz
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_inserter import DatabaseSourceInserter
from mirar.processors.database.database_selector import (
    CrossmatchSourceWithDatabase,
    DatabaseHistorySelector,
)
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.mask import (
    MaskAboveThreshold,
    MaskPixelsFromPath,
    MaskPixelsFromPathInverted,
    MaskPixelsFromWCS,
    WriteMaskedCoordsToFile,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.photometry.aperture_photometry import (
    ImageAperturePhotometry,
    SourceAperturePhotometry,
)
from mirar.processors.photometry.psf_photometry import (
    ImagePSFPhotometry,
    SourcePSFPhotometry,
)
from mirar.processors.reference import ProcessReference
from mirar.processors.sky import NightSkyMedianCalibrator
from mirar.processors.sources import CandidateNamer, SourceWriter, ZOGYSourceDetector
from mirar.processors.sources.source_table_builder import ForcedPhotometryCandidateTable
from mirar.processors.sources.source_table_modifier import CustomSourceModifier
from mirar.processors.sources.utils import RegionsWriter
from mirar.processors.utils import (
    HeaderAnnotator,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageSaver,
    ImageSelector,
)
from mirar.processors.utils.image_loader import LoadImageFromHeader
from mirar.processors.xmatch import XMatch
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [ImageLoader(input_sub_dir="raw", load_image=load_raw_wirc_image)]
# load_raw = [ImageLoader(input_sub_dir="firstpassstack",
# load_image=load_raw_wirc_image)]

log = [
    CSVLog(
        export_keys=[
            "OBJECT",
            "FILTER",
            "UTSHUT",
            "EXPTIME",
            "COADDS",
            OBSCLASS_KEY,
        ]
    )
]

masking = [MaskPixelsFromPath(mask_path=wirc_mask_path)]

dark_calibration = [ImageSelector(("exptime", "45.0")), DarkCalibrator()]

reduction = [
    ImageSaver(output_dir_name="darkcal"),
    HeaderAnnotator(input_keys=LATEST_SAVE_KEY, output_key=RAW_IMG_KEY),
    ImageDebatcher(),
    ImageSelector((OBSCLASS_KEY, "science")),
    # ImageSelector(("object", "ZTF18aavqmki")),
    ImageBatcher(split_key=["filter", "object"]),
    SkyFlatCalibrator(),
    NightSkyMedianCalibrator(),
    AutoAstrometry(catalog="tmc"),
    Sextractor(output_sub_dir="postprocess", **sextractor_astrometry_config),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
        cache=True,
    ),
    ImageSaver(output_dir_name="firstpass"),
    Swarp(swarp_config_path=swarp_sp_path, calculate_dims_in_swarp=True),
    ImageSaver(output_dir_name="firstpassstack"),
    # ImageSelector(("BASENAME", "image0125.fits_stack.fits")),
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
        copy_header_keys=FITS_MASK_KEY,
    ),
    LoadImageFromHeader(
        header_key=RAW_IMG_KEY,
        copy_header_keys=[SCAMP_HEADER_KEY, FITS_MASK_KEY],
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
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY),
    NightSkyMedianCalibrator(flat_mask_key=FITS_MASK_KEY),
    Sextractor(output_sub_dir="postprocess", **sextractor_astrometry_config),
    Swarp(swarp_config_path=swarp_sp_path, calculate_dims_in_swarp=True),
    Sextractor(output_sub_dir="final_sextractor", **sextractor_photometry_config),
    PhotCalibrator(
        ref_catalog_generator=wirc_photometric_catalog_generator,
        image_photometric_catalog_purifier=wirc_photometric_img_catalog_purifier,
        write_regions=True,
    ),
    ImageSaver(output_dir_name="final"),
]

reduce = log + masking + dark_calibration + reduction

reference = [
    ProcessReference(
        ref_image_generator=wirc_reference_image_generator,
        swarp_resampler=wirc_reference_image_resampler,
        sextractor=wirc_reference_sextractor,
        ref_psfex=wirc_reference_psfex,
    )
]

subtract = [
    Sextractor(**sextractor_reference_config, output_sub_dir="subtract", cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="subtract", norm_fits=True),
    ZOGYPrepare(output_sub_dir="subtract"),
    ZOGY(output_sub_dir="subtract", catalog_purifier=wirc_zogy_catalogs_purifier),
    ImageSaver(output_dir_name="diffs"),
]

image_photometry = [
    ImageAperturePhotometry(
        aper_diameters=[16],
        bkg_in_diameters=[25],
        bkg_out_diameters=[40],
        col_suffix_list=[""],
        phot_cutout_size=100,
        target_ra_key="TARGRA",
        target_dec_key="TARGDEC",
    ),
    Sextractor(**sextractor_reference_config, output_sub_dir="subtract", cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="photometry", norm_fits=True),
    ImagePSFPhotometry(target_ra_key="TARGRA", target_dec_key="TARGDEC"),
    ImageSaver(output_dir_name="photometry"),
]

export_candidates_from_header = [
    ForcedPhotometryCandidateTable(
        ra_header_key="TARGRA", dec_header_key="TARGDEC", name_header_key="TARGNAME"
    ),
]

candidate_photometry = [
    SourceAperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    SourcePSFPhotometry(),
]

detect_candidates = [
    ZOGYSourceDetector(
        output_sub_dir="subtract",
        **sextractor_candidate_config,
        copy_image_keywords=["PROGID", "PROGPI"],
    ),
]

process_candidates = [
    RegionsWriter(output_dir_name="candidates"),
    SourcePSFPhotometry(),
    SourceAperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    # SourceWriter(output_dir_name="candidates"),
    CustomSourceModifier(modifier_function=wirc_source_table_filter_annotator),
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
        db_output_columns=[CAND_NAME_KEY] + prv_candidate_cols,
    ),
    DatabaseSourceInserter(db_table=Candidate, duplicate_protocol="fail"),
    SourceWriter(output_dir_name="candidates"),
    # EdgeCandidatesMask(edge_boundary_size=100)
    # FilterCandidates(),
]

package_candidates = [
    IPACAvroExporter(
        base_name="WIRC",
        avro_schema_path=wirc_avro_schema_path,
    ),
    SendToFritz(
        base_name="WIRCTEST",
        group_ids=[1431],
        fritz_filter_id=74,
        instrument_id=5,
        stream_id=1005,
        update_thumbnails=True,
    ),
]

candidates = detect_candidates + process_candidates + package_candidates
imsub = reference + subtract + candidates
