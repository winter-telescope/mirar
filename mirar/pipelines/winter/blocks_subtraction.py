"""
Module with blocks for WINTER image subtraction and candidate processing
"""

from mirar.catalog.kowalski import PS1, TMASS
from mirar.paths import SOURCE_HISTORY_KEY, SOURCE_NAME_KEY, TARGET_KEY, ZP_KEY
from mirar.pipelines.winter.config import (
    prv_candidate_cols,
    psfex_path,
    sextractor_candidate_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
    winter_avro_schema_path,
    winter_fritz_config,
)
from mirar.pipelines.winter.generator import (
    mask_stamps_around_bright_stars,
    winter_candidate_annotator_filterer,
    winter_candidate_avro_fields_calculator,
    winter_candidate_quality_filterer,
    winter_history_deprecated_constraint,
    winter_imsub_catalog_purifier,
    winter_new_source_updater,
    winter_reference_generator,
    winter_reference_image_resampler_for_zogy,
    winter_reference_psf_phot_sextractor,
    winter_reference_psfex,
    winter_reference_sextractor,
    winter_skyportal_annotator,
    winter_source_entry_updater,
)
from mirar.pipelines.winter.models import (
    NAME_START,
    SOURCE_PREFIX,
    Candidate,
    Diff,
    Source,
)
from mirar.processors.astromatic import PSFex
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.avro import IPACAvroExporter
from mirar.processors.database.database_inserter import (
    DatabaseImageInserter,
    DatabaseSourceInserter,
)
from mirar.processors.database.database_selector import (
    SelectSourcesWithMetadata,
    SingleSpatialCrossmatchSource,
)
from mirar.processors.mask import MaskPixelsFromFunction
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import GetReferenceImage, ProcessReference
from mirar.processors.skyportal.skyportal_candidate import SkyportalCandidateUploader
from mirar.processors.sources import (
    CandidateNamer,
    CustomSourceTableModifier,
    SourceLoader,
    SourceWriter,
    ZOGYSourceDetector,
)
from mirar.processors.split import SUB_ID_KEY, SwarpImageSplitter
from mirar.processors.utils import (
    HeaderAnnotator,
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageSaver,
)
from mirar.processors.xmatch import XMatch
from mirar.processors.zogy.reference_aligner import AlignReference
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

# Reference building
refbuild = [
    GetReferenceImage(ref_image_generator=winter_reference_generator),
    ImageSaver(output_dir_name="stacked_ref"),
]

# Image subtraction
split_stack = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD", "STACKID"]),
    SwarpImageSplitter(swarp_config_path=swarp_config_path, n_x=2, n_y=1),
    ImageSaver(output_dir_name="split_stacks"),
]

imsub = [
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", TARGET_KEY, "SUBCOORD", "STACKID"]),
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

load_sources = [
    SourceLoader(input_dir_name="candidates"),
]

crossmatch_candidates = [
    XMatch(catalog=TMASS(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1(num_sources=3, search_radius_arcmin=0.5)),
    SourceWriter(output_dir_name="kowalski"),
    CustomSourceTableModifier(
        modifier_function=winter_candidate_avro_fields_calculator
    ),
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

name_candidates = (
    [
        # Check if the source is already in the source table
        SingleSpatialCrossmatchSource(
            db_table=Source,
            db_output_columns=["sourceid", SOURCE_NAME_KEY],
            crossmatch_radius_arcsec=2.0,
            ra_field_name="average_ra",
            dec_field_name="average_dec",
        ),
        # Assign names to the new sources
        CandidateNamer(
            db_table=Source,
            base_name=SOURCE_PREFIX,
            name_start=NAME_START,
            db_name_field=SOURCE_NAME_KEY,
        ),
        # Add the new sources to the source table
        CustomSourceTableModifier(modifier_function=winter_new_source_updater),
        DatabaseSourceInserter(
            db_table=Source,
            duplicate_protocol="ignore",
        ),
        # Get all candidates associated with source
    ]
    + select_history
    + [
        # Update average ra and dec for source
        CustomSourceTableModifier(modifier_function=winter_source_entry_updater),
        # Update sources in the source table
        DatabaseSourceInserter(
            db_table=Source,
            duplicate_protocol="replace",
        ),
        # Add candidates in the candidate table
        DatabaseSourceInserter(
            db_table=Candidate,
            duplicate_protocol="fail",
        ),
        SourceWriter(output_dir_name="preavro"),
    ]
)

avro_export = [
    IPACAvroExporter(
        topic_prefix="winter",
        base_name="WNTR",
        broadcast=False,
        avro_schema_path=winter_avro_schema_path,
    ),
    CustomSourceTableModifier(modifier_function=winter_candidate_quality_filterer),
    SourceWriter(output_dir_name="preskyportal"),
]

process_candidates = crossmatch_candidates + name_candidates + avro_export

load_skyportal = [SourceLoader(input_dir_name="preskyportal")]

send_to_skyportal = [
    CustomSourceTableModifier(modifier_function=winter_skyportal_annotator),
    SkyportalCandidateUploader(**winter_fritz_config),
]
