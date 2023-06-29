"""
Module containing WIRC-specific generator functions to
yield e.g catalog for astrometric calibrations
"""
import logging
import os

from mirar.catalog import Gaia2Mass
from mirar.data import Image
from mirar.pipelines.wirc.wirc_files import (
    psfex_path,
    sextractor_reference_config,
    wirc_file_dir,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.references.wirc import WIRCRef

logger = logging.getLogger(__name__)


def wirc_photometric_img_catalog_purifier(catalog, image):
    """
    Function to purify the photometric catalog

    :return: purified catalog
    """
    edge_width_pixels = 100
    fwhm_threshold_arcsec = 4.0

    x_lower_limit = edge_width_pixels
    x_upper_limit = image.get_data().shape[1] - edge_width_pixels
    y_lower_limit = edge_width_pixels
    y_upper_limit = image.get_data().shape[0] - edge_width_pixels

    clean_mask = (
        (catalog["FLAGS"] == 0)
        & (catalog["FWHM_WORLD"] < fwhm_threshold_arcsec / 3600.0)
        & (catalog["X_IMAGE"] > x_lower_limit)
        & (catalog["X_IMAGE"] < x_upper_limit)
        & (catalog["Y_IMAGE"] > y_lower_limit)
        & (catalog["Y_IMAGE"] < y_upper_limit)
    )

    return catalog[clean_mask]


def wirc_astrometric_catalog_generator(_) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for astrometry

    :return: catalogue
    """
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=10)


def wirc_photometric_catalog_generator(image: Image) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    filter_name = image["FILTER"]
    return Gaia2Mass(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=10,
        filter_name=filter_name,
        acceptable_h_ph_quals=["A"],
        acceptable_k_ph_quals=["A"],
    )


def wirc_reference_image_generator(
    image: Image,
    images_directory: str = os.environ.get("REF_IMG_DIR"),
) -> WIRCRef:
    """
    Function to match a new wirc image to a reference image directory

    :param image: image
    :param images_directory: ref image directory
    :return: wirc ref
    """
    object_name = image["OBJECT"]
    filter_name = image["FILTER"]
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory,
    )


def wirc_reference_image_resampler(**kwargs) -> Swarp:
    """Returns a SWarp resampler for WIRC"""
    return Swarp(
        swarp_config_path=wirc_file_dir.joinpath("config.swarp"), cache=True, **kwargs
    )


def wirc_reference_sextractor(output_sub_dir: str, gain: float) -> Sextractor:
    """Returns a Sextractor processor for WIRC reference images"""
    return Sextractor(
        **sextractor_reference_config,
        gain=gain,
        output_sub_dir=output_sub_dir,
        cache=True
    )


def wirc_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
    """Returns a PSFEx processor for WIRC"""
    return PSFex(
        config_path=psfex_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
    )
