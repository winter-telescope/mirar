"""
Module containing WIRC-specific generator functions to
yield e.g catalog for astrometric calibrations
"""
import logging
import os

import numpy as np
import pandas as pd
from astropy.table import Table

from mirar.catalog import Gaia2Mass
from mirar.data import Image
from mirar.paths import DIFF_IMG_KEY, FILTER_KEY, MAGLIM_KEY, REF_IMG_KEY, SCI_IMG_KEY
from mirar.pipelines.wirc.load_wirc_image import wirc_filter_dict
from mirar.pipelines.wirc.wirc_files import (
    psfex_path,
    sextractor_reference_config,
    wirc_file_dir,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.references.wirc import WIRCRef

logger = logging.getLogger(__name__)


def wirc_source_table_filter_annotator(src_df: pd.DataFrame) -> pd.DataFrame:
    """
    Function to remove bad candidates with None in sigmapsf, magpsf, magap, sigmagap,
    and update the source table with the keys required for the WIRC database
    :param src_df: Source dataframe
    :return: annotated dataframe
    """
    none_mask = (
        src_df.loc[:, "sigmapsf"].isnull()
        | src_df.loc[:, "magpsf"].isnull()
        | src_df.loc[:, "magap"].isnull()
        | src_df.loc[:, "sigmagap"].isnull()
    )

    mask = none_mask.values

    # Needing to do this because the dataframe is big-endian
    mask_inds = np.where(~mask)[0]
    src_df = pd.DataFrame([src_df.loc[x] for x in mask_inds]).reset_index(drop=True)

    src_df.loc[:, "sciimgname"] = src_df.loc[:, SCI_IMG_KEY]
    src_df.loc[:, "refimgname"] = src_df.loc[:, REF_IMG_KEY]
    src_df.loc[:, "diffimgname"] = src_df.loc[:, DIFF_IMG_KEY]
    src_df.loc[:, "fid"] = [wirc_filter_dict[x] for x in src_df.loc[:, FILTER_KEY]]
    src_df.loc[:, "programid"] = src_df.loc[:, "progid"]
    src_df.loc[:, "programpi"] = src_df.loc[:, "progpi"]
    src_df.loc[:, "diffmaglim"] = src_df.loc[:, MAGLIM_KEY]
    return src_df


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
        min_mag=10, max_mag=20, search_radius_arcmin=10, filter_name=filter_name
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
        swarp_config_path=wirc_file_dir.joinpath("config.swarp"),
        cache=True,
        subtract_bkg=True,
        **kwargs
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


def wirc_zogy_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table
) -> (Table, Table):
    """
    Function to purify the photometric catalog
    :param sci_catalog:
    :param ref_catalog:
    :return: sci_catalog, ref_catalog
    """
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["FLAGS"] == 0)
        & (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (ref_catalog["SNR_WIN"] < 1000)
    )
    return good_sci_sources, good_ref_sources
