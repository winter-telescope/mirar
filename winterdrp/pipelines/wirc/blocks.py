import os
from winterdrp.processors.dark import DarkCalibrator
from winterdrp.processors.flat import SkyFlatCalibrator
from winterdrp.processors.sky import NightSkyMedianCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path, sextractor_astrometry_config, scamp_fp_path, \
    swarp_sp_path
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
from winterdrp.processors.utils import ImageLoader
from winterdrp.processors.csvlog import CSVLog
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.generator import wirc_astrometric_catalog_generator, wirc_photometric_catalog_generator
from winterdrp.processors.alert_packets.avro_alert import AvroPacketMaker
from winterdrp.processors.send_to_fritz import SendToFritz
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher
from winterdrp.paths import core_fields, base_name_key
from winterdrp.processors.candidates.utils import RegionsWriter, DataframeWriter
from winterdrp.processors.photometry.psf_photometry import PSFPhotometry
from winterdrp.processors.photometry.aperture_photometry import AperturePhotometry
from winterdrp.catalog.kowalski import TMASS, PS1
from winterdrp.processors.xmatch import XMatch
from winterdrp.processors.database.database_importer import DatabaseHistoryImporter
from winterdrp.processors.database.postgres import get_colnames_from_schema
from winterdrp.processors.candidates.namer import CandidateNamer
from winterdrp.processors.database.database_exporter import DatabaseDataframeExporter
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.psfex import PSFex
from winterdrp.processors.reference import Reference
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare
from winterdrp.processors.candidates.candidate_detector import DetectCandidates
from winterdrp.pipelines.wirc.wirc_files import candidate_colnames, wirc_candidate_schema_path, \
    sextractor_reference_config, sextractor_candidate_config, psfex_path
from winterdrp.pipelines.wirc.generator import wirc_reference_image_generator, wirc_reference_sextractor, \
    wirc_reference_image_resampler, wirc_reference_psfex

load_raw = [ImageLoader(
        input_sub_dir="raw",
        load_image=load_raw_wirc_image
)]

reduce = [
    CSVLog(
        export_keys=["OBJECT", "FILTER", "UTSHUT", "EXPTIME", "COADDS", "OBSTYPE", "OBSCLASS"]
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
    Sextractor(
        output_sub_dir="postprocess",
        **sextractor_astrometry_config
    ),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
    ),
    Swarp(swarp_config_path=swarp_sp_path),
    Sextractor(
        output_sub_dir="final_sextractor",
        **sextractor_astrometry_config
    ),
    ImageSaver(output_dir_name="final"),
    PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator)
]

reference = [
    Reference(
        ref_image_generator=wirc_reference_image_generator,
        swarp_resampler=wirc_reference_image_resampler,
        sextractor=wirc_reference_sextractor,
        ref_psfex=wirc_reference_psfex
    )
]

subtract = [
    Sextractor(**sextractor_reference_config, output_sub_dir='subtract', cache=False),
    PSFex(config_path=psfex_path, output_sub_dir="subtract", norm_fits=True),
    ZOGYPrepare(output_sub_dir="subtract"),
    ZOGY(output_sub_dir="subtract"),
]

candidates = [
    DetectCandidates(output_sub_dir="subtract", **sextractor_candidate_config),
    RegionsWriter(output_dir_name='candidates'),
    PSFPhotometry(),
    AperturePhotometry(aper_diameters=[16, 70], cutout_size_aper_phot=100, bkg_in_diameters=[25, 90],
                       bkg_out_diameters=[40, 100], col_suffix_list=['', 'big']),
    DataframeWriter(output_dir_name='candidates'),
    XMatch(
        catalog=TMASS(),
        num_stars=3,
        search_radius_arcsec=30
    ),
    XMatch(
        catalog=PS1(),
        num_stars=3,
        search_radius_arcsec=30
    ),
    # History(),
    DataframeWriter(output_dir_name='kowalski'),
    DatabaseHistoryImporter(
        xmatch_radius_arcsec=2,
        time_field_name='jd',
        history_duration_days=500,
        db_name="wirc",
        db_user=os.environ.get('DB_USER'),
        db_pwd=os.environ.get('DB_PWD'),
        db_table='candidates',
        db_output_columns=candidate_colnames,
        schema_path=wirc_candidate_schema_path,
        q3c=False
    ),
    CandidateNamer(
        db_name='wirc',
        cand_table_name='candidates',
        base_name='WIRC',
        name_start='aaaaa',
        xmatch_radius_arcsec=2
    ),
    DatabaseDataframeExporter(
        db_name='wirc',
        db_table='candidates',
        schema_path=wirc_candidate_schema_path
    ),
    DataframeWriter(output_dir_name='dbop'),
    # EdgeCandidatesMask(edge_boundary_size=100)
    # FilterCandidates(),
    AvroPacketMaker(
        output_sub_dir="avro",
        base_name="WNTR",
        broadcast=False,
        save_local=True
    ),
    # SendToFritz(update_thumbnails = True)
]

imsub = reference + subtract + candidates
