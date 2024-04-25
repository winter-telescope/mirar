"""
Module containing functions to generate astrometric/photometric calibration catalogs
for WASP images
"""

# pylint: disable=duplicate-code

import logging

from mirar.catalog import BaseCatalog, Gaia2Mass
from mirar.catalog.vizier import PS1, SkyMapper
from mirar.catalog.vizier.sdss import SDSS, NotInSDSSError, in_sdss
from mirar.data.image_data import Image
from mirar.pipelines.git.config import (
    psfex_config_path,
    sextractor_photometry_config,
    swarp_config_path,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.processors.astromatic.sextractor.sextractor import SEXTRACTOR_HEADER_KEY
from mirar.references import BaseReferenceGenerator, PS1Ref, SDSSRef

logger = logging.getLogger(__name__)

WASP_SEARCH_RADIUS_ARCMIN = 30.0
WASP_PHOTOMETRIC_MAX_MAG = 22


def wasp_astrometric_catalog_generator(image: Image) -> Gaia2Mass:
    """
    Returns an astrometric catalog for summer,
    which is just a Gaia/2MASS one

    :param image: image to generate a catalog for
    :return: Gaia/2MASS catalog around image
    """
    temp_cat_path = image[SEXTRACTOR_HEADER_KEY]
    cat = Gaia2Mass(
        min_mag=10,
        max_mag=24,
        search_radius_arcmin=WASP_SEARCH_RADIUS_ARCMIN,
        trim=True,
        image_catalog_path=temp_cat_path,
        filter_name="j",
        acceptable_j_ph_quals=["A", "B", "C"],
    )
    return cat


def wasp_photometric_catalog_generator(image: Image) -> BaseCatalog:
    """
    Generate a photometric calibration catalog for SUMMER images

    For u band: SDSS if possible, otherwise Skymapper, otherwise fail
    For g/r/i/z: use PS1

    :param image: Image
    :return: catalog at image position
    """
    filter_name = image["FILTER"].replace("'", "")
    dec = image["OBJDEC"]

    if filter_name in ["u", "U"]:
        if in_sdss(image["OBJRA"], image["OBJDEC"]):
            return SDSS(
                min_mag=10,
                max_mag=WASP_PHOTOMETRIC_MAX_MAG,
                search_radius_arcmin=WASP_SEARCH_RADIUS_ARCMIN,
                filter_name=filter_name,
            )

        if dec < 0.0:
            return SkyMapper(
                min_mag=10,
                max_mag=WASP_PHOTOMETRIC_MAX_MAG,
                search_radius_arcmin=WASP_SEARCH_RADIUS_ARCMIN,
                filter_name=filter_name,
            )

        err = "U band image is in a field with no reference image."
        logger.error(err)
        raise NotInSDSSError(err)

    return PS1(
        min_mag=10,
        max_mag=WASP_PHOTOMETRIC_MAX_MAG,
        search_radius_arcmin=WASP_SEARCH_RADIUS_ARCMIN,
        filter_name=filter_name,
    )


def wasp_reference_image_generator(image: Image) -> BaseReferenceGenerator:
    """
    Get a reference image generator for a WASP image

    For u band: SDSS if possible, otherwise fail
    For g/r/i/z: use PS1

    :param image: image
    :return: Reference image generator
    """
    filter_name = image["FILTER"].replace("'", "")
    logger.debug(f"Filter is {filter_name}")

    if filter_name in ["u", "U"]:
        if in_sdss(image["CRVAL1"], image["CRVAL2"]):
            return SDSSRef(filter_name=filter_name)

        err = "U band image is in a field with no reference image."
        logger.error(err)
        raise NotInSDSSError(err)

    logger.debug("Will query reference image from PS1")
    return PS1Ref(filter_name=filter_name)


def wasp_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(
        swarp_config_path=swarp_config_path, cache=False, subtract_bkg=True, **kwargs
    )


def wasp_sdss_reference_cat_purifier(catalog, image: Image):
    """
    Purify SDSS catalog
    :param catalog:
    :param image:
    :return:
    """
    zero_point = image["ZP"]
    mags = catalog["MAG_AUTO"] + zero_point
    good_sources_mask = mags > 18
    return catalog[good_sources_mask]


def wasp_reference_sextractor(
    output_sub_dir: str, gain: float | None = None
) -> Sextractor:
    """
    Generates a sextractor processor for reference images

    :param output_sub_dir: output sui directory
    :param gain: gain of image
    :return: Sextractor processor
    """
    return Sextractor(
        gain=gain,
        output_sub_dir=output_sub_dir,
        cache=False,
        saturation=None,
        # catalog_purifier=git_sdss_reference_cat_purifier,
        **sextractor_photometry_config,
    )


def wasp_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
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


def wasp_zogy_catalogs_purifier(sci_catalog, ref_catalog):
    """
    Purify catalogs for ZOGY
    """
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
    )

    return good_sci_sources, good_ref_sources
