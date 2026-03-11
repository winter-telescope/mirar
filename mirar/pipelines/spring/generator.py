# from mirar.data import Image
import logging

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS

from mirar.catalog import PS1, Gaia2Mass
from mirar.data import Image, ImageBatch
from mirar.paths import FILTER_KEY, SEXTRACTOR_HEADER_KEY, get_output_dir
from mirar.pipelines.spring.config import (
    ref_psfex_path,
    sextractor_astrometry_config,
    sextractor_reference_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
)
from mirar.pipelines.spring.constants import spring_filters_map
from mirar.pipelines.spring.models import (
    RefComponent,
    RefQuery,
    RefStack,
    RefStacksTable,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.processors.base_catalog_xmatch_processor import (
    default_image_sextractor_catalog_purifier,
)
from mirar.references import PS1Ref
from mirar.references.wfcam.wfcam_query import WFAUQuery
from mirar.references.wfcam.wfcam_stack import WFCAMStackedRef
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


def spring_anet_sextractor_config_path_generator(_image: Image) -> str:
    """
    Generate the path to the ANET SExtractor configuration file for SPRING images
    Parameters
    ----------
    image:Image

    Returns
    -------
    sextractor_config_path:str

    """
    return sextractor_astrometry_config["config_path"]


def spring_photometric_catalog_generator(image: Image) -> Gaia2Mass | PS1:
    """
    Function to match SPRING image to GAIA/2MASS/PS1 for photometry
    Parameters
    ----------
    image:Image

    Returns
    ----------
    catalogue
    """
    filter_name = image[FILTER_KEY]
    search_radius_arcmin = (
        np.max([image["NAXIS1"], image["NAXIS2"]])
        * np.max([np.abs(image["CD1_1"]), np.abs(image["CD1_2"])])
        * 60.0
    ) / 2.0

    if filter_name in ["J", "H"]:
        return Gaia2Mass(
            min_mag=0,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name,
            snr_threshold=20,
            cache_catalog_locally=False,
        )

    if filter_name in ["Y"]:
        return PS1(
            min_mag=0,
            max_mag=20,
            search_radius_arcmin=search_radius_arcmin,
            filter_name=filter_name.lower(),
            cache_catalog_locally=False,
        )

    err = f"Filter name {filter_name} not recognized"
    logger.error(err)
    raise ValueError(err)


def spring_ref_photometric_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> tuple[Table, Table]:
    """
    Default function to purify the photometric image catalog
    """
    return default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog,
        ref_catalog=ref_catalog,
        image=image,
        edge_width_pixels=100,
        fwhm_threshold_arcsec=6,
    )


def spring_photcal_color_columns_generator(image: Image) -> tuple[list, list, tuple]:
    """
    Generates the columns for the Photometric calibration process.

    :param image: Image
    :return: (magnitude, error_magnitude)"""
    filter_name = image[FILTER_KEY]
    if filter_name in ["Y"]:
        return ["ymag", "zmag"], ["e_ymag", "e_zmag"], (0, 25)
    if filter_name == "J":
        return ["j_m", "h_m"], ["j_msigcom", "h_msigcom"], (0, 25)
    if filter_name == "H":
        return ["h_m", "k_m"], ["h_msigcom", "k_msigcom"], (0, 25)
    err = f"Filter name {filter_name} not recognized"
    logger.error(err)
    raise ValueError(err)


def spring_stackid_annotator(batch: ImageBatch) -> ImageBatch:
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


def spring_wfau_component_image_stacker(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """

    return Swarp(
        swarp_config_path=swarp_config_path,
        cache=False,
        include_scamp=False,
        combine=True,
        calculate_dims_in_swarp=True,
        subtract_bkg=True,
        center_type="ALL",
        **kwargs,
    )


def spring_reference_stack_annotator(stacked_image: Image, image: Image) -> Image:
    """
    Generates a stack id for SPRING reference images

    :param stacked_image: Image
    :param image: Image
    """
    stackid = (
        f"{str(image.header['FIELDID']).rjust(5, '0')}"
        f"{str(spring_filters_map[image.header['FILTER']])}"
    )
    stacked_image["STACKID"] = int(stackid)
    stacked_image["FIELDID"] = image.header["FIELDID"]
    return stacked_image


def spring_reference_generator(image: Image):
    """
    Gets a reference image generator for the spring data

    :param image: Image
    :return: Reference image generator
    """
    components_image_dir = get_output_dir(
        dir_root="components", sub_dir="spring/references"
    )
    components_image_dir.mkdir(parents=True, exist_ok=True)

    filtername = image["FILTER"]

    if filtername not in ["Y", "J", "H"]:
        raise ValueError(f"Filter {filtername} not recognized for SPRING")

    fieldid = int(image["FIELDID"])
    logger.debug(f"Fieldid: {fieldid}")

    cache_ref_stack = False
    if filtername in ["J", "H"]:
        wfcam_query = WFAUQuery(
            num_query_points=4,
            filter_name=filtername,
            use_db_for_component_queries=False,
            components_db_table=RefComponent,
            query_db_table=RefQuery,
            skip_online_query=False,
            component_image_subdir="spring/references/components",
        )

        return WFCAMStackedRef(
            filter_name=filtername,
            wfcam_query=wfcam_query,
            image_resampler_generator=spring_wfau_component_image_stacker,
            write_stacked_image=cache_ref_stack,
            write_stack_sub_dir="spring/references/ref_stacks",
            write_stack_to_db=cache_ref_stack,
            stacks_db_table=RefStack,
            component_image_sub_dir="components",
            references_base_subdir_name="spring/references",
            stack_image_annotator=spring_reference_stack_annotator,
        )

    assert filtername == "Y", f"Filter {filtername} not recognized for SPRING"

    logger.debug("Will query reference image from PS1")
    return PS1Ref(filter_name=filtername)


def spring_reference_image_resampler_for_zogy(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(
        swarp_config_path=swarp_config_path, cache=False, subtract_bkg=False, **kwargs
    )


def spring_reference_sextractor(
    output_sub_dir: str,
) -> Sextractor:
    """
    Returns a Sextractor processor for SPRING reference images

    :param output_sub_dir: str
    :return: Sextractor processor
    """
    return Sextractor(
        **sextractor_reference_config, output_sub_dir=output_sub_dir, cache=True
    )


def spring_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
    """
    Returns a PSFEx processor for SPRING

    :param output_sub_dir: str
    :param norm_fits: bool
    """
    return PSFex(
        config_path=ref_psfex_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
    )


def spring_reference_psf_phot_sextractor(output_sub_dir: str) -> Sextractor:
    """Returns a Sextractor processor for SPRING reference images

    :param output_sub_dir: str
    :return: Sextractor"""
    return Sextractor(
        **sextractor_reference_psf_phot_config,
        output_sub_dir=output_sub_dir,
        cache=False,
        use_psfex=True,
    )


def spring_imsub_catalog_purifier(sci_catalog: Table, ref_catalog: Table):
    """
    :param sci_catalog:
    :param ref_catalog:
    "return:
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
        & (ref_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (ref_catalog["SNR_WIN"] < 200)
        & (ref_catalog["FLUX_MAX"] < 30000)
    )

    return good_sci_sources, good_ref_sources


def mask_stamps_around_bright_stars(image: Image):
    """
    Masks the stamps around bright stars in the image
    :param image:
    :return: masked image"""

    catalog_path = image[SEXTRACTOR_HEADER_KEY]
    catalog = get_table_from_ldac(catalog_path)
    logger.debug(f"Using {catalog_path} as image catalog")
    bright_stars = catalog[(catalog["SNR_WIN"] > 1000)]
    logger.debug(f"Found {len(bright_stars)} bright stars in the image")
    stamp_half_size = 20
    bright_star_crds = SkyCoord(
        ra=bright_stars["ALPHAWIN_J2000"],
        dec=bright_stars["DELTAWIN_J2000"],
        unit="deg",
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
