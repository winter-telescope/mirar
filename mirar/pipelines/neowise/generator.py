import logging

from mirar.data import Image
from mirar.pipelines.neowise.config import (
    psfex_config_path,
    sextractor_reference_config,
    swarp_config_path,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.references.local import RefFromPath

logger = logging.getLogger(__name__)


def neowise_reference_image_generator(image: Image):
    """
    Gets a reference image generator for the wirc data

    :param image: Image
    :return: Reference image generator
    """
    filter_name = image["FILTER"]
    return RefFromPath(filter_name=filter_name, path=image["REF_PATH"])


def neowise_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(
        swarp_config_path=swarp_config_path, cache=True, subtract_bkg=True, **kwargs
    )


def neowise_reference_sextractor(
    output_sub_dir: str,
) -> Sextractor:
    """
    Generates a sextractor processor for reference images

    :param output_sub_dir: output sui directory
    :param gain: gain of image
    :return: Sextractor processor
    """
    return Sextractor(
        output_sub_dir=output_sub_dir,
        cache=True,
        saturation=10,
        # catalog_purifier=git_sdss_reference_cat_purifier,
        **sextractor_reference_config,
    )


def neowise_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
    """
    Generates a PSFex processor for reference images

    :param output_sub_dir: output sui directory
    :param norm_fits: boolean
    :return: Sextractor processor
    """
    return PSFex(
        config_path=psfex_config_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
    )


def neowise_zogy_catalogs_purifier(sci_catalog, ref_catalog):
    """
    Purify catalogs for ZOGY
    """
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 10.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 10.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
    )

    return good_sci_sources, good_ref_sources
