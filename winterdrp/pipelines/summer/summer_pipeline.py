import logging
import os
import astropy.io.fits
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.downloader.caltech import download_via_ssh
from astropy.time import Time
from winterdrp.catalog import Gaia2Mass, PS1, SDSS
from winterdrp.processors.database.database_exporter import DatabaseImageExporter
from winterdrp.processors.astromatic.sextractor.sextractor import sextractor_header_key
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp, PSFex
from astropy.io import fits
from winterdrp.pipelines.summer.summer_files import get_summer_schema_path, summer_mask_path, summer_weight_path, \
    sextractor_astrometry_config, sextractor_photometry_config, scamp_path, swarp_path
from winterdrp.pipelines.summer.summer_files.schema import summer_schema_dir
from winterdrp.processors.utils import ImageSaver
from winterdrp.processors.utils.image_loader import ImageLoader
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher
from winterdrp.processors.split import SplitImage, sub_id_key
from winterdrp.processors.utils import ImageSaver, HeaderAnnotator, ImageLoader, ImageSelector, ImageBatcher
from winterdrp.processors.utils.image_rejector import ImageRejector
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors import MaskPixels, BiasCalibrator, FlatCalibrator
from winterdrp.processors.csvlog import CSVLog
from winterdrp.paths import core_fields, base_name_key, latest_save_key, raw_img_key

# Image subtraction modules
from winterdrp.processors.reference import Reference
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare, default_summer_catalog_purifier
from winterdrp.references.ps1 import PS1Ref
from winterdrp.references.sdss import SDSSRef
from winterdrp.processors.candidates.candidate_detector import DetectCandidates
from winterdrp.processors.photometry.psf_photometry import PSFPhotometry
from winterdrp.processors.photometry.aperture_photometry import AperturePhotometry
from winterdrp.processors.candidates.utils import RegionsWriter, DataframeWriter
import pkg_resources

summer_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
summer_gain = 1.0
summer_pixel_scale = 0.466

logger = logging.getLogger(__name__)


def summer_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    temp_cat_path = header[sextractor_header_key]
    cat = Gaia2Mass(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=30,
        trim=True,
        image_catalog_path=temp_cat_path,
        filter_name='j'
    )
    return cat


def summer_photometric_catalog_generator(
        header: astropy.io.fits.Header
):
    filter_name = header['FILTERID']
    if filter_name == 'u':
        return SDSS(min_mag=10, max_mag=20, search_radius_arcmin=30, filter_name=filter_name)
    else:
        return PS1(min_mag=10, max_mag=20, search_radius_arcmin=30, filter_name=filter_name)


def load_raw_summer_image(
        path: str
) -> tuple[np.array, astropy.io.fits.Header]:
    with fits.open(path) as data:
        header = data[0].header
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "SCIENCE"]
        header['UTCTIME'] = header['UTCSHUT']
        header['TARGET'] = header['OBSTYPE'].lower()

        crd = SkyCoord(ra=data[0].header['RA'], dec=data[0].header['DEC'], unit=(u.deg, u.deg))
        header['RA'] = crd.ra.deg
        header['DEC'] = crd.dec.deg

        header['CRVAL1'] = header['RA']
        header['CRVAL2'] = header['DEC']

        tel_crd = SkyCoord(ra=data[0].header['TELRA'], dec=data[0].header['TELDEC'], unit=(u.deg, u.deg))
        header['TELRA'] = tel_crd.ra.deg
        header['TELDEC'] = tel_crd.dec.deg
        header['BZERO'] = 0

        header[latest_save_key] = path
        header[raw_img_key] = path

        data[0].data = data[0].data * 1.0

        if 'other' in header['FILTERID']:
            header['FILTERID'] = 'r'

        header["CALSTEPS"] = ""

        base_name = os.path.basename(path)
        header[base_name_key] = base_name
        header["EXPID"] = int("".join(base_name.split("_")[1:3])[2:])
        # header["EXPID"] = str(header["NIGHT"]) + str(header["OBSHISTID"])

        pipeline_version = pkg_resources.require("winterdrp")[0].version
        pipeline_version_padded_str = "".join([x.rjust(2, "0") for x in pipeline_version.split(".")])
        header["PROCID"] = int(str(header["EXPID"])+str(pipeline_version_padded_str))
        header.append(('GAIN', summer_gain, 'Gain in electrons / ADU'), end=True)

        header['OBSDATE'] = int(header['UTC'].split('_')[0])

        obstime = Time(header['UTCISO'], format='iso')
        t0 = Time('2018-01-01', format='iso')
        header['NIGHT'] = np.floor((obstime - t0).jd).astype(int)

        header['EXPMJD'] = header['OBSMJD']

        for key in ["PROGID", "OBSID"]:
            if key not in header.keys():
                # logger.warning(f"No {key} found in header of {path}")
                header[key] = 0

        if "SUBPROG" not in header.keys():
            # logger.warning(f"No SUBPROG found in header of {path}")
            header['SUBPROG'] = 'none'

        header['FILTER'] = header['FILTERID']
        try:
            header['SHUTOPEN'] = Time(header['SHUTOPEN'], format='iso').jd
        except (KeyError, ValueError):
            logger.warning(f"Error parsing 'SHUTOPEN' of {path}: ({header['SHUTOPEN']})")

        try:
            header['SHUTCLSD'] = Time(header['SHUTCLSD'], format='iso').jd
        except ValueError:
            logger.warning(f"Error parsing 'SHUTCLSD' of {path}: ({header['SHUTCLSD']})")

        header['PROCFLAG'] = 0
        sunmoon_keywords = [
            'MOONRA', 'MOONDEC', 'MOONILLF', 'MOONPHAS', 'MOONALT', 'SUNAZ', 'SUNALT',
        ]
        for key in sunmoon_keywords:
            val = 0
            if key in header.keys():
                if header[key] not in ['']:
                    val = header[key]
            header[key] = val

        itid_dict = {
            'SCIENCE': 1,
            'BIAS': 2,
            'FLAT': 2,
            'DARK': 2,
            'FOCUS': 3,
            'POINTING': 4,
            'OTHER': 5
        }

        if not header['OBSTYPE'] in itid_dict.keys():
            header['ITID'] = 5
        else:
            header['ITID'] = itid_dict[header['OBSTYPE']]

        if header['FIELDID'] == 'radec':
            header['FIELDID'] = 999999999

        if header['ITID'] != 1:
            header['FIELDID'] = -99

        if 'COADDS' not in header.keys():
            header['COADDS'] = 1

        if header['PROGID'] == '' or header['PROGID'] == 'WINTER':
            header['PROGID'] = 2
        crds = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=(u.deg, u.deg))
        header['RA'] = crds.ra.deg
        header['DEC'] = crds.dec.deg

        data[0].header = header
    return data[0].data, data[0].header


def load_proc_summer_image(
        path: str
) -> tuple[np.array, astropy.io.fits.Header]:
    img = fits.open(path)
    data = img[0].data
    header = img[0].header
    if 'ZP' not in header.keys():
        header['ZP'] = header['ZP_AUTO']
        header['ZP_std'] = header['ZP_AUTO_std']

    header['CENTRA'] = header['CRVAL1']
    header['CENTDEC'] = header['CRVAL2']
    pipeline_version = pkg_resources.require("winterdrp")[0].version
    pipeline_version_padded_str = "".join([x.rjust(2, "0") for x in pipeline_version.split(".")])
    header['DIFFID'] = int(str(header["EXPID"])+str(pipeline_version_padded_str))
    data[data == 0] = np.nan
    # logger.info(header['CRVAL2'])
    return data, header


def summer_reference_image_generator(
        header: fits.header,
):
    filter_name = header['FILTER']
    logger.info(f'Filter is {filter_name}')
    if filter_name == 'u':
        logger.info(f'Will query reference image from SDSS')
        return SDSSRef(filter_name=filter_name)
    else:
        logger.info(f'Will query reference image from PS1')
        return PS1Ref(filter_name=filter_name)


def summer_reference_image_resampler(pixscale,
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
    return Swarp(swarp_config_path='winterdrp/pipelines/summer/summer_imsub_files/config/config.swarp',
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


def summer_reference_sextractor(output_sub_dir, gain):
    return Sextractor(config_path='winterdrp/pipelines/summer/summer_imsub_files/config/photomCat.sex',
                      parameter_path='winterdrp/pipelines/summer/summer_imsub_files/config/photom.param',
                      filter_path='winterdrp/pipelines/summer/summer_imsub_files/config/default.conv',
                      starnnw_path='winterdrp/pipelines/summer/summer_imsub_files/config/default.nnw',
                      gain=gain,
                      output_sub_dir=output_sub_dir,
                      cache=True
                      )


def summer_reference_psfex(output_sub_dir, norm_fits):
    return PSFex(config_path='winterdrp/pipelines/summer/summer_imsub_files/config/photom.psfex',
                 output_sub_dir=output_sub_dir,
                 norm_fits=norm_fits,
                 cache=True
                 )


pipeline_name = "summer"


class SummerPipeline(Pipeline):
    name = pipeline_name

    all_pipeline_configurations = {
        None: [
            ImageLoader(
                load_image=load_raw_summer_image
            ),
            CSVLog(
                export_keys=[
                                "UTC", 'FIELDID', "FILTERID", "EXPTIME", "OBSTYPE", "RA", "DEC", "TARGTYPE",
                                base_name_key
                            ] + core_fields
            ),
            ImageRejector(("OBSTYPE", "FOCUS"), ("FILTER", "?")),
            DatabaseImageExporter(
                db_name=pipeline_name,
                db_table="exposures",
                schema_path=get_summer_schema_path("exposures"),
                full_setup=True,
                schema_dir=summer_schema_dir
            ),
            MaskPixels(mask_path=summer_mask_path),
            DatabaseImageExporter(
                db_name=pipeline_name,
                db_table="raw",
                schema_path=get_summer_schema_path("raw")
            ),
            BiasCalibrator(),
            ImageBatcher(split_key="filter"),
            FlatCalibrator(),
            ImageBatcher(base_name_key),
            ImageSelector(("OBSTYPE", "SCIENCE")),
            ImageSelector((base_name_key, "SUMMER_20220816_042349_Camera0.fits")),
            # ImageSelector((base_name_key, "SUMMER_20220402_214324_Camera0.fits")),
            AutoAstrometry(pa=0, inv=True, pixel_scale=summer_pixel_scale),
            # ImageLoader(input_sub_dir='autoastrometry',
            #             load_image=load_raw_summer_image),
            # ImageSelector((base_name_key, "SUMMER_20220816_041023_Camera0.fits")),
            Sextractor(
                output_sub_dir="testb",
                weight_image=summer_weight_path,
                checkimage_name=None,
                checkimage_type=None,
                **sextractor_astrometry_config
            ),
            Scamp(
                ref_catalog_generator=summer_astrometric_catalog_generator,
                scamp_config_path=scamp_path,
            ),
            Swarp(swarp_config_path=swarp_path, imgpixsize=2400),
            # ImageSaver(output_dir_name="photprocess"),
            Sextractor(output_sub_dir="photprocess",
                       checkimage_name='NONE',
                       checkimage_type='NONE',
                       **sextractor_photometry_config),
            # ImageSaver(output_dir_name="processed"),
            PhotCalibrator(ref_catalog_generator=summer_photometric_catalog_generator),
            ImageSaver(output_dir_name="processed", additional_headers=['PROCIMG'], write_mask=True),
            DatabaseImageExporter(
                db_name=pipeline_name,
                db_table="proc",
                schema_path=get_summer_schema_path("proc")
            )
        ],
        'imsub': [
            ImageLoader(
                input_sub_dir='processed',
                load_image=load_proc_summer_image
            ),
            ImageBatcher(split_key=base_name_key),
            ImageSelector(('OBSTYPE', 'SCIENCE')),
            # ImageSelector(('FILTER', ['u'])),
            ImageSelector((base_name_key, "SUMMER_20220816_042349_Camera0.resamp.fits")),
            Reference(ref_image_generator=summer_reference_image_generator,
                      ref_psfex=summer_reference_psfex,
                      ref_sextractor=summer_reference_sextractor,
                      ref_swarp_resampler=summer_reference_image_resampler),
            Sextractor(config_path='winterdrp/pipelines/summer/summer_imsub_files/config/photomCat.sex',
                       parameter_path='winterdrp/pipelines/summer/summer_imsub_files/config/photom.param',
                       filter_path='winterdrp/pipelines/summer/summer_imsub_files/config/default.conv',
                       starnnw_path='winterdrp/pipelines/summer/summer_imsub_files/config/default.nnw',
                       output_sub_dir='subtract',
                       cache=False,
                       write_regions_file=True),
            PSFex(config_path='winterdrp/pipelines/summer/summer_imsub_files/config/photom.psfex',
                  output_sub_dir="subtract",
                  norm_fits=True),
            ImageSaver(output_dir_name='ref'),
            ZOGYPrepare(output_sub_dir="subtract", sci_zp_header_key='ZP_AUTO', catalog_purifier=default_summer_catalog_purifier),
            ZOGY(output_sub_dir="subtract"),
            DatabaseImageExporter(
                db_name=pipeline_name,
                db_table="diff",
                schema_path=get_summer_schema_path("diff"),

            ),
            DetectCandidates(output_sub_dir="subtract",
                             cand_det_sextractor_config='winterdrp/pipelines/summer/summer_imsub_files/config/photomCat.sex',
                             cand_det_sextractor_nnw='winterdrp/pipelines/summer/summer_imsub_files/config/default.nnw',
                             cand_det_sextractor_filter='winterdrp/pipelines/summer/summer_imsub_files/config/default.conv',
                             cand_det_sextractor_params='winterdrp/pipelines/summer/summer_imsub_files/config/Scorr.param'),
            RegionsWriter(output_dir_name='candidates'),
            PSFPhotometry(),
            AperturePhotometry(aper_diameters=[8, 40], cutout_size_aper_phot=100, bkg_in_diameters=[25, 90],
                               bkg_out_diameters=[40, 100], col_suffix_list=['', 'big']),
            DataframeWriter(output_dir_name='candidates'),
        ]
    }

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        download_via_ssh(
            server="jagati.caltech.edu",
            base_dir="/data/viraj/winter_data/commissioning/raw/",
            night=night,
            pipeline=pipeline_name
        )
