"""
Module with generators for WINTER pipeline
"""
import logging
import os
from typing import Type

import numpy as np

from mirar.catalog import Gaia2Mass
from mirar.data import Image
from mirar.paths import get_output_dir
from mirar.pipelines.winter.models import RefComponents, RefQueries, RefStacks
from mirar.pipelines.wirc.wirc_files import sextractor_astrometry_config
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.references.local import RefFromPath
from mirar.references.ukirt import UKIRTRef

logger = logging.getLogger(__name__)
winter_dir = os.path.dirname(__file__)
astromatic_config_dir = os.path.join(winter_dir, "config/")
swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")


def winter_reference_generator(image: Image, db_table: Type[BaseDB] = RefStacks):
    """
    Generates a reference image for the winter data
    Args:
        db_table: Database table to search for existing image
        image: Image

    Returns:

    """
    components_image_dir = get_output_dir(
        dir_root="components", sub_dir="winter/references"
    )
    if not components_image_dir.exists():
        components_image_dir.mkdir(parents=True)

    filtername = image["FILTER"]
    # TODO if in_ukirt and in_vista, different processing
    fieldid = int(image["FIELDID"])
    subdetid = int(image["SUBDETID"])
    logger.debug(f"Fieldid: {fieldid}, subdetid: {subdetid}")
    db_results = db_table.sql_model().select_query(
        select_keys=["savepath"],
        compare_keys=["fieldid", "subdetid"],
        compare_values=[fieldid, subdetid],
        comparators=["__eq__", "__eq__"],
    )

    if len(db_results) > 0:
        savepaths = [x[0] for x in db_results]
        if os.path.exists(savepaths[0]):
            logger.info(f"Found reference image in database: {savepaths[0]}")
            return RefFromPath(path=savepaths[0], filter_name=filtername)

    return UKIRTRef(
        filter_name=filtername,
        swarp_resampler=winter_reference_image_resampler,
        sextractor_generator=ref_sextractor,
        phot_calibrator_generator=winter_reference_phot_calibrator,
        num_query_points=9,
        query_table=RefQueries,
        components_table=RefComponents,
        write_to_db=True,
        write_db_table=RefStacks,
        component_image_dir=components_image_dir.as_posix(),
        night_sub_dir="winter/references",
    )


def winter_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    logger.info(kwargs)
    return Swarp(
        swarp_config_path=swarp_config_path, subtract_bkg=True, cache=False, **kwargs
    )


def winter_photometric_catalog_generator(image: Image) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    filter_name = image["FILTER"]
    search_radius_arcmin = (
        np.max([image["NAXIS1"], image["NAXIS2"]]) * np.abs(image["CD1_1"]) * 60
    )
    return Gaia2Mass(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=search_radius_arcmin,
        filter_name=filter_name,
        snr_threshold=20,
    )


def winter_reference_phot_calibrator(image: Image, **kwargs) -> PhotCalibrator:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    x_lower_limit = 0
    y_lower_limit = 0
    x_upper_limit = image.header["NAXIS1"]
    y_upper_limit = image.header["NAXIS2"]

    return PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        x_lower_limit=x_lower_limit,
        x_upper_limit=x_upper_limit,
        y_lower_limit=y_lower_limit,
        y_upper_limit=y_upper_limit,
        write_regions=True,
        **kwargs,
    )


def ref_sextractor(image: Image):
    """
    Generates a sextractor instance for reference images to get photometry
    Args:
        image:

    Returns:

    """
    logger.debug(image)
    return Sextractor(
        output_sub_dir="phot",
        **sextractor_astrometry_config,
        write_regions_bool=True,
        cache=False,
    )


def ref_phot_calibrator(image: Image):
    """
    Generates a photcalibrator instance for reference images to get photometry
    Args:
        image:

    Returns:

    """
    logger.debug(image)
    return PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        write_regions=True,
        fwhm_threshold_arcsec=3,
    )
