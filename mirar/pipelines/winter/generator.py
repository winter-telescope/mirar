"""
Module with generators for WINTER pipeline
"""
import logging
import os

import numpy as np
from astropy.table import Table

from mirar.catalog import Gaia2Mass
from mirar.catalog.vizier import PS1
from mirar.data import Image
from mirar.paths import SATURATE_KEY, get_output_dir
from mirar.pipelines.winter.config import (
    psfex_path,
    sextractor_reference_config,
    swarp_config_path,
)
from mirar.pipelines.winter.constants import winter_filters_map
from mirar.pipelines.winter.fourier_bkg_model import subtract_fourier_background_model
from mirar.pipelines.winter.models import RefComponents, RefQueries, RefStacks
from mirar.pipelines.wirc.wirc_files import sextractor_astrometry_config
from mirar.processors.astromatic import PSFex
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.base_catalog_xmatch_processor import (
    default_image_sextractor_catalog_purifier,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.split import SUB_ID_KEY
from mirar.references.local import RefFromPath
from mirar.references.ukirt import UKIRTRef

logger = logging.getLogger(__name__)


def winter_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    logger.debug(kwargs)
    return Swarp(swarp_config_path=swarp_config_path, cache=False, **kwargs)


def winter_reference_sextractor(output_sub_dir: str, gain: float) -> Sextractor:
    """Returns a Sextractor processor for WINTER reference images"""
    return Sextractor(
        **sextractor_reference_config,
        gain=gain,
        output_sub_dir=output_sub_dir,
        cache=True,
    )


def winter_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
    """Returns a PSFEx processor for WINTER"""
    return PSFex(
        config_path=psfex_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
    )


def winter_astrostat_catalog_purifier(catalog: Table, image: Image) -> Table:
    """
    Default function to purify the photometric image catalog
    """

    return default_image_sextractor_catalog_purifier(
        catalog, image, edge_width_pixels=0, fwhm_threshold_arcsec=20.0
    )


def winter_photometric_catalog_generator(image: Image) -> Gaia2Mass | PS1:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    filter_name = image["FILTER"]
    search_radius_arcmin = (
        np.max([image["NAXIS1"], image["NAXIS2"]])
        * np.max([np.abs(image["CD1_1"]), np.abs(image["CD1_2"])])
        * 60
    ) / 2.0

    if filter_name in ["J", "H"]:
        return Gaia2Mass(
            min_mag=10,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name,
            snr_threshold=20,
        )

    if filter_name in ["Y"]:
        return PS1(
            min_mag=10,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name.lower(),
        )

    err = f"Filter {filter_name} not recognised"
    logger.error(err)
    raise ValueError(err)


def winter_ref_photometric_img_catalog_purifier(catalog: Table, image: Image) -> Table:
    """
    Default function to purify the photometric image catalog
    """

    return default_image_sextractor_catalog_purifier(
        catalog, image, edge_width_pixels=100, fwhm_threshold_arcsec=4.0
    )


def winter_reference_phot_calibrator(_: Image, **kwargs) -> PhotCalibrator:
    """
    Generates a resampler for reference images

    :param _: image
    :param kwargs: kwargs
    :return: Swarp processor
    """

    return PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        write_regions=True,
        image_photometric_catalog_purifier=winter_ref_photometric_img_catalog_purifier,
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


def winter_astrometric_ref_catalog_generator(_) -> Gaia2Mass:
    """
    Function to generate a reference catalog for WINTER astrometry

    :return: catalogue
    """
    return Gaia2Mass(min_mag=7, max_mag=20, search_radius_arcmin=20)


def winter_astrometry_sextractor_catalog_purifier(catalog: Table, _) -> Table:
    """
    Function to purify the Sextractor catalog for WINTER astrometry
    """
    clean_catalog = catalog[
        (catalog["FLAGS"] == 0) & (catalog["FWHM_IMAGE"] > 0) & (catalog["SNR_WIN"] > 0)
    ]
    return clean_catalog


def winter_stackid_annotator(image: Image) -> Image:
    """
    Generates a stack id for WINTER images

    :param image: Image
    :return: stack id
    """
    first_rawid = np.min([int(x) for x in image["RAWID"].split(",")])
    image["STACKID"] = int(first_rawid)
    return image


def winter_reference_stackid_generator(image: Image) -> int:
    """
    Generates a stack id for WINTER reference images
    """
    stackid = (
        f"{str(image.header['FIELDID']).rjust(5, '0')}"
        f"{str(image.header[SUB_ID_KEY]).rjust(2, '0')}"
        f"{str(winter_filters_map[image.header['FILTER']])}"
    )
    return int(stackid)


def winter_reference_generator(image: Image):
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
    subdetid = int(image[SUB_ID_KEY])
    logger.debug(f"Fieldid: {fieldid}, subdetid: {subdetid}")
    # db_results = db_table.sql_model().select_query(
    #     select_keys=["savepath"],
    #     compare_keys=["fieldid", SUB_ID_KEY.lower()],
    #     compare_values=[fieldid, subdetid],
    #     comparators=["__eq__", "__eq__"],
    # )
    #
    # if len(db_results) > 0:
    #     savepaths = [x[0] for x in db_results]
    #     if os.path.exists(savepaths[0]):
    #         logger.debug(f"Found reference image in database: {savepaths[0]}")
    #         return RefFromPath(path=savepaths[0], filter_name=filtername)

    return UKIRTRef(
        filter_name=filtername,
        swarp_resampler=winter_reference_image_resampler,
        sextractor_generator=ref_sextractor,
        phot_calibrator_generator=winter_reference_phot_calibrator,
        num_query_points=9,
        query_table=RefQueries,
        components_table=RefComponents,
        write_to_db=False,
        write_db_table=RefStacks,
        component_image_dir=components_image_dir.as_posix(),
        night_sub_dir="winter/references",
        skip_online_query=False,
        stack_id_generator=winter_reference_stackid_generator,
    )


def winter_refbuild_reference_generator(image: Image):
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
    stack_image_dir = get_output_dir(
        dir_root="stacked_ref", sub_dir="winter/references"
    )
    if not components_image_dir.exists():
        components_image_dir.mkdir(parents=True)

    filtername = image["FILTER"]
    # TODO if in_ukirt and in_vista, different processing
    fieldid = int(image["FIELDID"])
    subdetid = int(image[SUB_ID_KEY])
    logger.debug(f"Fieldid: {fieldid}, subdetid: {subdetid}")

    stack_image_path = stack_image_dir / f"field{fieldid}_subdet{subdetid}_ref.fits"

    if os.path.exists(stack_image_path):
        logger.debug(f"Found reference image in database: {stack_image_path}")
        return RefFromPath(path=stack_image_path, filter_name=filtername)

    return UKIRTRef(
        filter_name=filtername,
        swarp_resampler=winter_reference_image_resampler,
        sextractor_generator=ref_sextractor,
        phot_calibrator_generator=winter_reference_phot_calibrator,
        num_query_points=9,
        query_table=RefQueries,
        components_table=RefComponents,
        write_to_db=False,
        write_db_table=RefStacks,
        component_image_dir=components_image_dir.as_posix(),
        night_sub_dir="winter/references",
        skip_online_query=False,
        stack_id_generator=winter_reference_stackid_generator,
    )


def winter_fourier_filtered_image_generator(image: Image) -> Image:
    """
    Generates a fourier filtered image for the winter data
    """
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

    return image
