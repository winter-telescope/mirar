import logging
import os

import astropy
from astropy.io import fits

from winterdrp.catalog import Gaia2Mass
from winterdrp.pipelines.wirc.wirc_files import (
    psfex_path,
    sextractor_reference_config,
    wirc_file_dir,
)
from winterdrp.processors.astromatic import PSFex, Sextractor, Swarp
from winterdrp.references.wirc import WIRCRef

logger = logging.getLogger(__name__)


def wirc_astrometric_catalog_generator(header: astropy.io.fits.Header) -> Gaia2Mass:
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30)


def wirc_photometric_catalog_generator(header: astropy.io.fits.Header) -> Gaia2Mass:
    filter_name = header["FILTER"]
    return Gaia2Mass(
        min_mag=10, max_mag=20, search_radius_arcmin=30, filter_name=filter_name
    )


def wirc_reference_image_generator(
    header: fits.header,
    images_directory: str = os.environ.get("REF_IMG_DIR"),
) -> WIRCRef:
    object_name = header["OBJECT"]
    filter_name = header["FILTER"]
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory,
    )


def wirc_reference_image_resampler(**kwargs) -> Swarp:
    return Swarp(
        swarp_config_path=wirc_file_dir.joinpath("config.swarp"), cache=True, **kwargs
    )


def wirc_reference_sextractor(output_sub_dir: str, gain: float) -> Sextractor:
    return Sextractor(
        **sextractor_reference_config,
        gain=gain,
        output_sub_dir=output_sub_dir,
        cache=True
    )


def wirc_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
    return PSFex(
        config_path=psfex_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
        cache=True,
    )
