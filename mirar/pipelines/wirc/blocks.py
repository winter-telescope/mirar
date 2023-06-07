"""
Module containing standard processing blocks for WIRC
"""
from mirar.catalog.kowalski import PS1, TMASS
from mirar.pipelines.wirc.generator import (
    wirc_astrometric_catalog_generator,
    wirc_photometric_catalog_generator,
    wirc_reference_image_generator,
    wirc_reference_image_resampler,
    wirc_reference_psfex,
    wirc_reference_sextractor,
)
from mirar.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from mirar.pipelines.wirc.wirc_files import (
    candidate_colnames,
    psfex_path,
    scamp_fp_path,
    sextractor_astrometry_config,
    sextractor_candidate_config,
    sextractor_reference_config,
    swarp_sp_path,
    wirc_candidate_schema_path,
    wirc_mask_path,
)
from mirar.processors.alert_packets.avro_alert import AvroPacketMaker
from mirar.processors.astromatic import Scamp, Sextractor, Swarp
from mirar.processors.astromatic.psfex import PSFex
from mirar.processors.autoastrometry import AutoAstrometry
from mirar.processors.candidates.candidate_detector import DetectCandidates
from mirar.processors.candidates.candidate_extractor import (
    ForcedPhotometryCandidateTable,
)
from mirar.processors.candidates.namer import CandidateNamer
from mirar.processors.candidates.utils import DataframeWriter, RegionsWriter
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.database.database_exporter import DatabaseDataframeExporter
from mirar.processors.database.database_importer import DatabaseHistoryImporter
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.mask import MaskPixels
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.photometry.aperture_photometry import (
    CandidateAperturePhotometry,
    ImageAperturePhotometry,
)
from mirar.processors.photometry.psf_photometry import (
    CandidatePSFPhotometry,
    ImagePSFPhotometry,
)
from mirar.processors.reference import ProcessReference
from mirar.processors.sky import NightSkyMedianCalibrator
from mirar.processors.utils import ImageLoader, ImageSaver
from mirar.processors.utils.image_selector import (
    ImageBatcher,
    ImageDebatcher,
    ImageSelector,
)
from mirar.processors.xmatch import XMatch
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [ImageLoader(input_sub_dir="raw", load_image=load_raw_wirc_image)]

reduce = [
    CSVLog(
        export_keys=[
            "OBJECT",
            "FILTER",
            "UTSHUT",
            "EXPTIME",
            "COADDS",
            "OBSTYPE",
            "OBSCLASS",
        ]
    ),
    MaskPixels(mask_path=wirc_mask_path),
    ImageSelector(("exptime", "45.0")),
    DarkCalibrator(),
    ImageDebatcher(),
    ImageSelector(("obsclass", "science")),
    ImageBatcher(split_key="filter"),
    SkyFlatCalibrator(),
    NightSkyMedianCalibrator(),
    AutoAstrometry(catalog="tmc"),
    Sextractor(output_sub_dir="postprocess", **sextractor_astrometry_config),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
    ),
    Swarp(swarp_config_path=swarp_sp_path),
    Sextractor(output_sub_dir="final_sextractor", **sextractor_astrometry_config),
    PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator),
    ImageSaver(output_dir_name="final"),
]

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
    ZOGY(output_sub_dir="subtract"),
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
    CandidateAperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    CandidatePSFPhotometry(),
]

detect_candidates = [
    DetectCandidates(output_sub_dir="subtract", **sextractor_candidate_config)
]

process_candidates = [
    RegionsWriter(output_dir_name="candidates"),
    CandidatePSFPhotometry(),
    CandidateAperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    DataframeWriter(output_dir_name="candidates"),
    XMatch(catalog=TMASS(num_sources=3, search_radius_arcmin=0.5)),
    XMatch(catalog=PS1(num_sources=3, search_radius_arcmin=0.5)),
    DataframeWriter(output_dir_name="kowalski"),
    DatabaseHistoryImporter(
        crossmatch_radius_arcsec=2.0,
        time_field_name="jd",
        history_duration_days=500.0,
        db_name="wirc",
        db_table="candidates",
        db_output_columns=candidate_colnames,
        schema_path=wirc_candidate_schema_path,
        q3c_bool=False,
    ),
    CandidateNamer(
        db_name="wirc",
        db_table="candidates",
        base_name="WIRC",
        name_start="aaaaa",
        xmatch_radius_arcsec=2,
        schema_path=wirc_candidate_schema_path,
    ),
    DatabaseDataframeExporter(
        db_name="wirc", db_table="candidates", schema_path=wirc_candidate_schema_path
    ),
    DataframeWriter(output_dir_name="dbop"),
    # EdgeCandidatesMask(edge_boundary_size=100)
    # FilterCandidates(),
]

package_candidates = [
    AvroPacketMaker(
        output_sub_dir="avro", base_name="WNTR", broadcast=False, save_local=True
    )
    # SendToFritz(update_thumbnails = True)
]

candidates = detect_candidates + process_candidates + package_candidates
imsub = reference + subtract + candidates
