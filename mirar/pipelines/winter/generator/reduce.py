"""
Module with generators for WINTER pipeline
"""

import logging
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS

from mirar.data import Image, ImageBatch
from mirar.database.constraints import DBQueryConstraints
from mirar.errors.exceptions import ProcessorError
from mirar.paths import (
    FILTER_KEY,
    OBSCLASS_KEY,
    REF_CAT_PATH_KEY,
    SATURATE_KEY,
    get_output_dir,
)
from mirar.pipelines.winter.config import sextractor_anet_config
from mirar.pipelines.winter.fourier_bkg_model import subtract_fourier_background_model
from mirar.processors.split import SUB_ID_KEY
from mirar.processors.utils.image_selector import select_from_images
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


class ReductionQualityError(ProcessorError):
    """Error raised when the quality of the reduction is too poor"""


def winter_stackid_annotator(batch: ImageBatch) -> ImageBatch:
    """
    Generates a stack id for WINTER images as the minimum of the RAWID of the
    images for which the stack was requested.

    :param batch: ImageBatch
    :return: ImageBatch with stackid added to the header
    """
    first_rawid = np.min([int(image["RAWID"]) for image in batch])
    for image in batch:
        image["STACKID"] = int(first_rawid)
    return batch


winter_history_deprecated_constraint = DBQueryConstraints(
    columns="deprecated", accepted_values="False", comparison_types="="
)


def winter_fourier_filtered_image_generator(batch: ImageBatch) -> ImageBatch:
    """
    Generates a fourier filtered image for the winter data
    """
    new_batch = []
    for image in batch:
        # First, set the nans in the raw_data to the median value
        raw_data = image.get_data()
        replace_value = np.nanmedian(raw_data)  # 0.0

        mask = image.get_mask()  # 0 is masked, 1 is unmasked

        raw_data[~mask] = replace_value

        filtered_data, sky_model = subtract_fourier_background_model(raw_data)

        # mask the data back
        filtered_data[~mask] = np.nan

        image.set_data(filtered_data)

        # Update the header
        image.header["MEDCOUNT"] = np.nanmedian(filtered_data)
        image.header[SATURATE_KEY] -= np.nanmedian(sky_model)
        new_batch.append(image)
    new_batch = ImageBatch(new_batch)
    return new_batch


def select_winter_sky_flat_images(images: ImageBatch) -> ImageBatch:
    """
    Selects the flat for the winter data, get the top 250 images sorted by median counts
    """
    flat_images, medcounts = [], []
    for image in images:
        image["MEDCOUNT"] = np.nanmedian(image.get_data())
        if image["MEDCOUNT"] > 2000:
            flat_images.append(image)
            medcounts.append(image["MEDCOUNT"])

    sort_inds = np.argsort(medcounts)
    flat_images = [flat_images[i] for i in sort_inds[::-1]][:250]
    flat_images = ImageBatch(flat_images)

    if len(flat_images) == 0:
        # To enable current version of realtime processing?
        logger.warning(
            "No good flat images found, using all images in batch. "
            f"The filter and subdetid of the first image in batch is "
            f"{images[0]['FILTER']} and {images[0][SUB_ID_KEY]}"
        )
        flat_images = select_from_images(
            images, key=OBSCLASS_KEY, target_values="science"
        )
    return flat_images


def select_winter_flat_images(images: ImageBatch) -> ImageBatch:
    """
    Selects the flat for the winter data, get the top 250 images sorted by median counts
    """
    flat_images = select_from_images(images, key=OBSCLASS_KEY, target_values="flat")
    return flat_images


def select_winter_dome_flats_images(images: ImageBatch) -> ImageBatch:
    """
    Selects the flat for the winter data, get the top 250 images sorted by median counts
    """
    flat_images = select_from_images(images, key=OBSCLASS_KEY, target_values="flat")
    flat_images = ImageBatch(
        [image for image in flat_images if image["MEDCOUNT"] > 20000.0]
    )
    return flat_images


def winter_master_flat_path_generator(images: ImageBatch) -> Path:
    """
    Generates a master flat path for the winter data

    :param images:
    :return: Path to master flat
    """
    filters_list = [image[FILTER_KEY] for image in images]
    image_filter = np.unique(filters_list)
    assert len(image_filter) == 1, "More than one filter in batch"
    image_filter = image_filter[0]
    subdetid_list = [image[SUB_ID_KEY] for image in images]
    subdetid = np.unique(subdetid_list)
    assert len(subdetid) == 1, "More than one subdetid in batch"
    subdetid = subdetid[0]

    master_flat_dir = get_output_dir(dir_root="winter/master_calibrations/masterflats")

    master_flat_path = master_flat_dir / f"master_flat_{image_filter}_{subdetid}.fits"
    return master_flat_path


def winter_anet_sextractor_config_path_generator(image: Image) -> str:
    """
    Generates the sextractor config file path for the winter image
    """
    if image["BOARD_ID"] in [1, 5, 6]:
        return sextractor_anet_config["config_path_boardid_1_5_6"]

    return sextractor_anet_config["config_path_boardid_2_3_4"]


def winter_imsub_catalog_purifier(sci_catalog: Table, ref_catalog: Table):
    """

    :param sci_catalog:
    :param ref_catalog:
    :return:
    """
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["FLAGS_MODEL"] == 0)
        & (sci_catalog["SNR_WIN"] > 10)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 200)
        & (sci_catalog["FLUX_MAX"] < 30000)
    )

    good_ref_sources = (
        (ref_catalog["FLAGS"] == 0)
        & (ref_catalog["FLAGS_MODEL"] == 0)
        & (ref_catalog["SNR_WIN"] > 10)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (ref_catalog["SNR_WIN"] < 1000)
        & (ref_catalog["FLUX_MAX"] < 30000)
    )

    return good_sci_sources, good_ref_sources


def mask_stamps_around_bright_stars(image: Image):
    """
    Masks the stamps around bright stars in the image
    :param image:
    :return: masked image
    """

    catalog_path = image[REF_CAT_PATH_KEY]
    catalog = get_table_from_ldac(catalog_path)
    bright_stars = catalog[(catalog["magnitude"] < 11)]
    logger.debug(f"Found {len(bright_stars)} bright stars in the image")
    stamp_half_size = 20
    bright_star_crds = SkyCoord(
        ra=bright_stars["ra"], dec=bright_stars["dec"], unit="deg"
    )
    wcs = WCS(image.get_header())
    bright_star_pix_x, bright_star_pix_y = wcs.all_world2pix(
        bright_star_crds.ra, bright_star_crds.dec, 1
    )
    logger.debug(f"Masking stamps around {len(bright_star_pix_x)} bright stars")
    mask = np.zeros_like(image.get_data(), dtype=bool)

    for x, y in zip(bright_star_pix_x, bright_star_pix_y):
        mask[
            int(y) - stamp_half_size : int(y) + stamp_half_size,
            int(x) - stamp_half_size : int(x) + stamp_half_size,
        ] = True

    return mask


def winter_boardid_6_demasker(images: ImageBatch) -> ImageBatch:
    """
    Demasks images from board 6 by replacing the bad channel pixels with the median of
    the unmasked pixels. This is required because swarp does not handle masked pixels
    distributed across the image well, producing a fully masked image.
    :param images: ImageBatch
    :return: ImageBatch
    """
    for image in images:
        boardid = image.header["BOARD_ID"]
        if boardid == 6:
            img_data = image.get_data()
            img_data[0::2, 0::4] = np.nanmedian(img_data)
            image.set_data(img_data)

    return images
