from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.references.wirc import WIRCRef
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.psfex import PSFex
from winterdrp.processors.reference import Reference
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare
from winterdrp.processors.candidates.candidate_detector import DetectCandidates
from astropy.io import fits, ascii
import os
import logging

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
from winterdrp.pipelines.wirc.wirc_pipeline import load_raw_wirc_image
from winterdrp.processors.database.database_importer import DatabaseHistoryImporter
from winterdrp.processors.database.postgres import get_colnames_from_schema
from winterdrp.processors.candidates.namer import CandidateNamer
from winterdrp.processors.database.database_exporter import DatabaseDataframeExporter
logger = logging.getLogger(__name__)


def wirc_reference_image_generator(
        header: fits.header,
        images_directory: str = os.environ.get('REF_IMG_DIR'),
):
    object_name = header['OBJECT']
    filter_name = header['FILTER']
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory
    )


def wirc_reference_image_resampler(pixscale,
                                   x_imgpixsize,
                                   y_imgpixsize,
                                   center_ra,
                                   center_dec,
                                   propogate_headerlist,
                                   temp_output_sub_dir,
                                   night_sub_dir,
                                   include_scamp,
                                   combine,
                                   gain,
                                   subtract_bkg):
    logger.debug(f'Night sub dir is {night_sub_dir}')
    return Swarp(swarp_config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/config.swarp',
                 pixscale=pixscale,
                 x_imgpixsize=x_imgpixsize,
                 y_imgpixsize=y_imgpixsize,
                 center_ra=center_ra,
                 center_dec=center_dec,
                 propogate_headerlist=propogate_headerlist,
                 temp_output_sub_dir=temp_output_sub_dir,
                 night_sub_dir=night_sub_dir,
                 include_scamp=include_scamp,
                 combine=combine,
                 gain=gain,
                 cache=True,
                 subtract_bkg=subtract_bkg
                 )


def wirc_reference_sextractor(output_sub_dir, gain):
    return Sextractor(config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photomCat.sex',
                      parameter_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photom.param',
                      filter_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.conv',
                      starnnw_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.nnw',
                      gain=gain,
                      output_sub_dir=output_sub_dir,
                      cache=True
                      )


def wirc_reference_psfex(output_sub_dir, norm_fits):
    return PSFex(config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photom.psfex',
                 output_sub_dir=output_sub_dir,
                 norm_fits=norm_fits,
                 cache=True
                 )


def detect_candidates_sextractor():
    pass


class WircImsubPipeline(Pipeline):
    name = "wirc_imsub"

    header_keys = [
        "UTSHUT",
        'OBJECT',
        "FILTER",
        "EXPTIME"
    ]
    batch_split_keys = ["UTSHUT"]

    candidates_db_columns = get_colnames_from_schema('winterdrp/pipelines/wirc_imsub/wirc_imsub_files/schema'
                                                     '/candidates.sql')
    all_pipeline_configurations = {
        None: [
            ImageLoader(
                input_sub_dir="raw",
                load_image=load_raw_wirc_image
            ),
            # ImageBatcher(split_key='UTSHUT'),
            ImageSelector((base_name_key, "ZTF21aagppzg_J_stack_1_20210702.fits")),
            Reference(
                ref_image_generator=wirc_reference_image_generator,
                ref_swarp_resampler=wirc_reference_image_resampler,
                ref_sextractor=wirc_reference_sextractor,
                ref_psfex=wirc_reference_psfex
            ),
            # Swarp(),
            Sextractor(config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photomCat.sex',
                       parameter_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photom.param',
                       filter_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.conv',
                       starnnw_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.nnw',
                       output_sub_dir='subtract',
                       cache=False),
            PSFex(config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photom.psfex',
                  output_sub_dir="subtract",
                  norm_fits=True),
            ZOGYPrepare(output_sub_dir="subtract"),
            ZOGY(output_sub_dir="subtract"),
            DetectCandidates(output_sub_dir="subtract",
                             cand_det_sextractor_config='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/photomCat.sex',
                             cand_det_sextractor_nnw='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.nnw',
                             cand_det_sextractor_filter='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/default.conv',
                             cand_det_sextractor_params='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/Scorr.param'),
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
                db_output_columns=candidates_db_columns,
                schema_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/schema/candidates.sql',
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
                schema_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/schema/candidates.sql'
            ),
            DataframeWriter(output_dir_name='dbop')
            # EdgeCandidatesMask(edge_boundary_size=100)
            # FilterCandidates(),
            # AvroPacketMaker(output_sub_dir="avro",
                            # base_name="WNTR",
                            # broadcast=False,
                            # save_local=False),
            SendToFritz()
        ]
    }
