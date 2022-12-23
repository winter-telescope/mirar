import astropy
import logging
from astropy.io import fits

from winterdrp.paths import sextractor_header_key
from winterdrp.references.ps1 import PS1Ref
from winterdrp.references.sdss import SDSSRef
from winterdrp.catalog import Gaia2Mass, PS1, SkyMapper
from winterdrp.catalog.sdss import SDSS, in_sdss, NotInSDSSError
from winterdrp.processors.astromatic.sextractor.sextractor import sextractor_header_key
from winterdrp.processors.astromatic import Sextractor, Swarp, PSFex
from winterdrp.pipelines.sedmv2.config import swarp_config_path, sextractor_photometry_config, psfex_config_path

logger = logging.getLogger(__name__)


def sedmv2_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    temp_cat_path = header[sextractor_header_key]
    cat = Gaia2Mass(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=7.5,
        trim=True,
        image_catalog_path=temp_cat_path,
        filter_name='j'
    )
    return cat


def sedmv2_photometric_catalog_generator(
        header: astropy.io.fits.Header
):
    filter_name = header['FILTERID']
    dec = header['DEC']
    if filter_name in ['u', "U"]:
        if in_sdss(header["RA"], header["DEC"]):
            return SDSS(min_mag=10, max_mag=20, search_radius_arcmin=7.5, filter_name=filter_name)
        elif dec < 0.:
            return SkyMapper(min_mag=10, max_mag=20, search_radius_arcmin=7.5, filter_name=filter_name)
        else:
            err = "U band image is in a field with no reference image."
            logger.error(err)
            raise NotInSDSSError(err)
    else:
        return PS1(min_mag=10, max_mag=20, search_radius_arcmin=7.5, filter_name=filter_name)


def sedmv2_reference_image_generator(
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


def sedmv2_reference_image_resampler(
        **kwargs
) -> Swarp:
    return Swarp(
        swarp_config_path=swarp_config_path,
        cache=True,
        **kwargs
        )


def sedmv2_reference_sextractor(output_sub_dir, gain):
    return Sextractor(
        gain=gain,
        output_sub_dir=output_sub_dir,
        cache=True,
        **sextractor_photometry_config
        )


def sedmv2_reference_psfex(output_sub_dir, norm_fits):
    return PSFex(
        config_path=psfex_config_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
        cache=True
        )