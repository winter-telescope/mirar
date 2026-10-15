"""
Module containing functions to generate astrometric/photometric calibration catalogs
for LMI images
"""

# pylint: disable=duplicate-code

import logging

from astropy.table import Table

from mirar.catalog import BaseCatalog, GaiaVizier
from mirar.catalog.vizier import PS1, SkyMapper
from mirar.catalog.vizier.sdss import SDSS, NotInSDSSError, in_sdss
from mirar.data.image_data import Image
from mirar.paths import BASE_NAME_KEY
from mirar.pipelines.git.config import (
    psfex_config_path,
    sextractor_photometry_config,
    swarp_config_path,
)
from mirar.pipelines.lmi.config.constants import LMI_WIDTH_DEG
from mirar.pipelines.lmi.generator.sources import lmi_skyportal_formatter
from mirar.pipelines.lmi.generator.stacks import label_stack_id
from mirar.pipelines.lmi.generator.target import annotate_target_coordinates
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.processors.base_catalog_xmatch_processor import (
    default_image_sextractor_catalog_purifier,
)
from mirar.references import BaseReferenceGenerator, PS1Ref, SDSSRef

logger = logging.getLogger(__name__)

LMI_SEARCH_RADIUS_ARCMIN = LMI_WIDTH_DEG * 60.0 / (2.0**0.5)  # Triangle diagonal
LMI_PHOTOMETRIC_MAX_MAG = 22


def lmi_astrometric_catalog_generator(_) -> GaiaVizier:
    """
    Returns an astrometric catalog for LMI

    :param _: Image
    :return: Gaia catalog around image
    """
    cat = GaiaVizier(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=LMI_SEARCH_RADIUS_ARCMIN,
        filter_name="G",
    )
    return cat


def lmi_photometric_catalog_generator(image: Image) -> BaseCatalog:
    """
    Generate a photometric calibration catalog for LMI images

    For u band: SDSS if possible, otherwise Skymapper, otherwise fail
    For g/r/i/z: use PS1

    :param image: Image
    :return: catalog at image position
    """
    filter_name = image["FILTER"].replace("'", "")

    ra, dec = image["CRVAL1"], image["CRVAL2"]

    if filter_name in ["u", "U"]:
        if in_sdss(ra, dec):
            return SDSS(
                min_mag=8,
                max_mag=LMI_PHOTOMETRIC_MAX_MAG,
                search_radius_arcmin=LMI_SEARCH_RADIUS_ARCMIN,
                filter_name=filter_name,
            )

        if dec < 0.0:
            return SkyMapper(
                min_mag=8,
                max_mag=LMI_PHOTOMETRIC_MAX_MAG,
                search_radius_arcmin=LMI_SEARCH_RADIUS_ARCMIN,
                filter_name=filter_name,
            )

        err = "U band image is in a field with no reference image."
        logger.error(err)
        raise NotInSDSSError(err)

    return PS1(
        min_mag=10,
        max_mag=LMI_PHOTOMETRIC_MAX_MAG,
        search_radius_arcmin=LMI_SEARCH_RADIUS_ARCMIN,
        filter_name=filter_name,
    )


def lmi_reference_image_generator(image: Image) -> BaseReferenceGenerator:
    """
    Get a reference image generator for a LMI image

    For u band: SDSS if possible, otherwise fail
    For g/r/i/z: use PS1

    :param image: image
    :return: Reference image generator
    """
    filter_name = image["FILTER"].replace("'", "")
    logger.debug(f"Filter is {filter_name}")

    if filter_name in ["u", "U"]:
        return SDSSRef(filter_name=filter_name)
    logger.debug("Will query reference image from PS1")
    return PS1Ref(filter_name=filter_name)


def lmi_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(
        swarp_config_path=swarp_config_path, cache=False, subtract_bkg=True, **kwargs
    )


def lmi_sdss_reference_cat_purifier(catalog, image: Image):
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


def lmi_reference_sextractor(
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


def lmi_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
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


def lmi_zogy_catalogs_purifier(sci_catalog, ref_catalog):
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


# Catalogue matches are dominated by galaxies, and MAG_AUTO grows with an
# extended source by an amount which differs source to source
LMI_CLASS_STAR_THRESHOLD = 0.8
LMI_ZP_SNR_THRESHOLD = 20.0

# Seeing bad enough to smear every star into the extended class would otherwise
# leave too few sources to fit a zeropoint at all
LMI_MIN_ZP_POINT_SOURCES = 10


def lmi_photometric_catalog_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> tuple[Table, Table]:
    """
    Keep only well-measured point sources for the photometric calibration

    The default purifier cuts on flags, seeing and the image edge, which leaves a
    catalogue in which galaxies outnumber stars, and those scatter the zeropoint
    several times more than the stars do.

    :param sci_catalog: Catalog of sources detected in the image
    :param ref_catalog: Catalog of reference sources
    :param image: Image being calibrated
    :return: Purified image and reference catalogs
    """
    sci_catalog, ref_catalog = default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog, ref_catalog=ref_catalog, image=image
    )

    point_sources = (sci_catalog["CLASS_STAR"] > LMI_CLASS_STAR_THRESHOLD) & (
        sci_catalog["SNR_WIN"] > LMI_ZP_SNR_THRESHOLD
    )

    if point_sources.sum() < LMI_MIN_ZP_POINT_SOURCES:
        logger.warning(
            f"Found only {point_sources.sum()} point sources in "
            f"{image[BASE_NAME_KEY]}, so the zeropoint is fit to all "
            f"{len(sci_catalog)} sources and will be less precise"
        )
        return sci_catalog, ref_catalog

    return sci_catalog[point_sources], ref_catalog
