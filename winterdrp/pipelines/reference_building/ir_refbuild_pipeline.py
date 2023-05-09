"""
Module with pipeline for building reference images in the IR from WFAU
"""
import logging
import os

from winterdrp.catalog import Gaia2Mass
from winterdrp.data import Image
from winterdrp.paths import get_output_dir
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.pipelines.reference_building.db_models import RefComponents, RefStacks
from winterdrp.pipelines.wirc.wirc_files import sextractor_astrometry_config
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.reference import GetReferenceImage
from winterdrp.processors.utils import ImageDebatcher, ImageSaver
from winterdrp.references.ukirt import UKIRTRef

logger = logging.getLogger(__name__)

refbuild_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(refbuild_dir, "config/")
swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")


def winter_reference_generator(image: Image):
    """
    Generates a reference image for the winter data
    Args:
        image: Image

    Returns:

    """
    filtername = image["FILTER"]
    # TODO check if exists in DB
    # TODO if in_ukirt and in_vista, different processing
    components_image_dir = get_output_dir(
        dir_root="components", sub_dir="ir_reference_building" "/references"
    )
    if not components_image_dir.exists():
        components_image_dir.mkdir(parents=True)

    return UKIRTRef(
        filter_name=filtername,
        swarp_resampler=winter_reference_image_resampler,
        sextractor_generator=ref_sextractor,
        phot_calibrator_generator=winter_reference_phot_calibrator,
        num_query_points=16,
        components_table=RefComponents,
        write_to_db=True,
        write_db_table=RefStacks,
        component_image_dir=components_image_dir.as_posix(),
        night_sub_dir="ir_reference_building/references",
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
    x_upper_limit = image.header["NAXIS2"]
    y_upper_limit = image.header["NAXIS1"]

    return PhotCalibrator(
        ref_catalog_generator=wirc_photometric_catalog_generator,
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
        ref_catalog_generator=wirc_photometric_catalog_generator,
        write_regions=True,
        fwhm_threshold_arcsec=3,
    )


class IRRefBuildPipeline(Pipeline):
    """
    Pipeline for building reference images in the IR from WFAU
    """

    name = "ir_reference_building"

    refbuild = [
        ImageDebatcher(),
        GetReferenceImage(
            ref_image_generator=winter_reference_generator,
        ),
        ImageSaver(output_dir_name="stacked_ref"),
    ]

    all_pipeline_configurations = {"default": refbuild}

    gain = 1.0
    non_linear_level = 65535

    @staticmethod
    def _load_raw_image(path: str):
        pass

    @staticmethod
    def download_raw_images_for_night(night: str):
        pass
