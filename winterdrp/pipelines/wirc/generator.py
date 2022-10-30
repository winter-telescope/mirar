import astropy
import os
import logging
from astropy.io import fits
from winterdrp.catalog import Gaia2Mass
from winterdrp.references.wirc import WIRCRef
from winterdrp.processors.astromatic import Sextractor, Swarp, PSFex

logger = logging.getLogger(__name__)


def wirc_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30)


def wirc_photometric_catalog_generator(
        header: astropy.io.fits.Header
):
    filter_name = header['FILTER']
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30, filter_name=filter_name)


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


def wirc_reference_image_resampler(**kwargs):
    return Swarp(swarp_config_path='winterdrp/pipelines/wirc_imsub/wirc_imsub_files/config/config.swarp',
                 cache=True, **kwargs
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