"""
Module containing functions to generate astrometric/photometric calibration catalogs
for SUMMER
"""
import logging

from winterdrp.catalog import PS1, BaseCatalog, Gaia2Mass, SkyMapper
from winterdrp.catalog.sdss import SDSS, NotInSDSSError, in_sdss
from winterdrp.data.image_data import Image
from winterdrp.pipelines.summer.config import (
    psfex_config_path,
    sextractor_photometry_config,
    swarp_config_path,
)
from winterdrp.processors.astromatic import PSFex, Sextractor, Swarp
from winterdrp.processors.astromatic.sextractor.sextractor import SEXTRACTOR_HEADER_KEY
from winterdrp.references.ps1 import PS1Ref
from winterdrp.references.sdss import SDSSRef

logger = logging.getLogger(__name__)


def summer_astrometric_catalog_generator(image: Image) -> BaseCatalog:
    temp_cat_path = image[SEXTRACTOR_HEADER_KEY]
    cat = Gaia2Mass(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=7.5,
        trim=True,
        image_catalog_path=temp_cat_path,
        filter_name="j",
    )
    return cat


def summer_photometric_catalog_generator(image: Image) -> BaseCatalog:
    filter_name = image["FILTERID"]
    dec = image["DEC"]

    if filter_name in ["u", "U"]:
        if in_sdss(image["RA"], image["DEC"]):
            return SDSS(
                min_mag=10,
                max_mag=20,
                search_radius_arcmin=7.5,
                filter_name=filter_name,
            )

        if dec < 0.0:
            return SkyMapper(
                min_mag=10,
                max_mag=20,
                search_radius_arcmin=7.5,
                filter_name=filter_name,
            )

        err = "U band image is in a field with no reference image."
        logger.error(err)
        raise NotInSDSSError(err)

    return PS1(
        min_mag=10, max_mag=20, search_radius_arcmin=7.5, filter_name=filter_name
    )


def summer_reference_image_generator(
    image: Image,
):
    filter_name = image["FILTER"]
    logger.info(f"Filter is {filter_name}")

    if filter_name in ["u", "U"]:
        logger.info("Will query reference image from SDSS")

        return SDSSRef(filter_name=filter_name)

    logger.info("Will query reference image from PS1")
    return PS1Ref(filter_name=filter_name)


def summer_reference_image_resampler(**kwargs) -> Swarp:
    return Swarp(swarp_config_path=swarp_config_path, cache=True, **kwargs)


def summer_reference_sextractor(output_sub_dir, gain):
    return Sextractor(
        gain=gain,
        output_sub_dir=output_sub_dir,
        cache=True,
        **sextractor_photometry_config,
    )


def summer_reference_psfex(output_sub_dir, norm_fits):
    return PSFex(
        config_path=psfex_config_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
    )
