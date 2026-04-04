# from mirar.data import Image
import logging

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS

from mirar.catalog import PS1, Gaia2Mass
from mirar.data import Image, ImageBatch
from mirar.data.source_data import SourceBatch
from mirar.errors.exceptions import ProcessorError
from mirar.paths import (
    FILTER_KEY,
    MAGLIM_KEY,
    SEXTRACTOR_HEADER_KEY,
    TIME_KEY,
    ZP_KEY,
    ZP_STD_KEY,
    get_output_dir,
)
from mirar.pipelines.mirage.config import (
    MIRAGE_GAIN,
    ref_psfex_path,
    sextractor_astrometry_config,
    sextractor_reference_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
)
from mirar.pipelines.mirage.constants import mirage_filters_map
from mirar.pipelines.mirage.models import (
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


class NoGoodCandidatesError(ProcessorError):
    """Error raised when no candidates pass quality cuts"""


def mirage_anet_sextractor_config_path_generator(_image: Image) -> str:
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


def mirage_photometric_catalog_generator(image: Image) -> Gaia2Mass | PS1:
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


def mirage_ref_photometric_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table, image: Image
) -> tuple[Table, Table]:
    """
    Default function to purify the photometric image catalog
    """
    return default_image_sextractor_catalog_purifier(
        sci_catalog=sci_catalog,
        ref_catalog=ref_catalog,
        image=image,
        edge_width_pixels=0,
        fwhm_threshold_arcsec=6,
    )


def mirage_photcal_color_columns_generator(image: Image) -> tuple[list, list, tuple]:
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


def mirage_stackid_annotator(batch: ImageBatch) -> ImageBatch:
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


def mirage_stack_gain_modifier(batch: ImageBatch) -> ImageBatch:
    """
    Ad-Hoc function for SPRING Images to ensure that the GAIN keyword is
    written in properly instead of being left blank

    :param batch: ImageBatch
    :return: ImageBatch with GAIN added to the header
    """
    coadds = batch[0]["COADDS"]
    effective_gain = coadds * MIRAGE_GAIN
    for image in batch:
        if image["GAIN"] <= MIRAGE_GAIN:
            image["GAIN"] = effective_gain
    return batch


def mirage_wfau_component_image_stacker(**kwargs) -> Swarp:
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


def mirage_reference_stack_annotator(stacked_image: Image, image: Image) -> Image:
    """
    Generates a stack id for SPRING reference images

    :param stacked_image: Image
    :param image: Image
    """
    stackid = (
        f"{str(image.header['FIELDID']).rjust(5, '0')}"
        f"{str(mirage_filters_map[image.header['FILTER']])}"
    )
    stacked_image["STACKID"] = int(stackid)
    stacked_image["FIELDID"] = image.header["FIELDID"]
    return stacked_image


def mirage_reference_generator(image: Image):
    """
    Gets a reference image generator for the mirage data

    :param image: Image
    :return: Reference image generator
    """
    components_image_dir = get_output_dir(
        dir_root="components", sub_dir="mirage/references"
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
            component_image_subdir="mirage/references/components",
        )

        return WFCAMStackedRef(
            filter_name=filtername,
            wfcam_query=wfcam_query,
            image_resampler_generator=mirage_wfau_component_image_stacker,
            write_stacked_image=cache_ref_stack,
            write_stack_sub_dir="mirage/references/ref_stacks",
            write_stack_to_db=cache_ref_stack,
            stacks_db_table=RefStack,
            component_image_sub_dir="components",
            references_base_subdir_name="mirage/references",
            stack_image_annotator=mirage_reference_stack_annotator,
        )

    assert filtername == "Y", f"Filter {filtername} not recognized for SPRING"

    logger.debug("Will query reference image from PS1")
    return PS1Ref(filter_name=filtername)


def mirage_reference_image_resampler_for_zogy(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(
        swarp_config_path=swarp_config_path, cache=False, subtract_bkg=False, **kwargs
    )


def mirage_reference_sextractor(
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


def mirage_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
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


def mirage_reference_psf_phot_sextractor(output_sub_dir: str) -> Sextractor:
    """Returns a Sextractor processor for SPRING reference images

    :param output_sub_dir: str
    :return: Sextractor"""
    return Sextractor(
        **sextractor_reference_psf_phot_config,
        output_sub_dir=output_sub_dir,
        cache=False,
        use_psfex=True,
    )


def mirage_imsub_catalog_purifier(sci_catalog: Table, ref_catalog: Table):
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
        & (ref_catalog["SNR_WIN"] < 1000)
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


def mirage_candidate_annotator_filterer(source_batch: SourceBatch) -> SourceBatch:
    """
    Function to perform basic filtering to weed out bad candidates with None
    magnitudes, to be added.
    :param source_batch: Source batch
    :return: updated batch
    """

    new_batch = []

    for source in source_batch:
        src_df = source.get_data()

        bad_sources_mask = (
            src_df["sigmapsf"].isnull()
            | src_df["magpsf"].isnull()
            | src_df["magap"].isnull()
            | src_df["sigmagap"].isnull()
            | (src_df["fwhm"] <= 0)
            | (src_df["scorr"] < 0)
        )

        mask = bad_sources_mask.values

        # Needing to do this because the dataframe is big-endian
        mask_inds = np.where(~mask)[0]
        filtered_df = pd.DataFrame([src_df.loc[x] for x in mask_inds]).reset_index(
            drop=True
        )

        if len(filtered_df) == 0:
            filtered_df = pd.DataFrame(columns=src_df.columns)

        # Pipeline (db) specific keywords
        source["magzpsci"] = source[ZP_KEY]
        source["magzpsciunc"] = source[ZP_STD_KEY]
        source["diffmaglim"] = source[MAGLIM_KEY]
        source["programpi"] = source["PROGPI"]
        source["programid"] = source["PROGID"]
        source["field"] = source["FIELDID"]
        tstr = source[TIME_KEY]
        new_tstr = (
            tstr[:4]
            + "-"
            + tstr[4:6]
            + "-"
            + tstr[6:8]
            + "T"
            + tstr[9:11]
            + ":"
            + tstr[11:13]
            + ":"
            + tstr[13:]
        )
        source["jd"] = Time(new_tstr).jd

        source.set_data(filtered_df)
        if len(filtered_df) > 0:
            new_batch.append(source)

    if len(new_batch) == 0:
        raise NoGoodCandidatesError

    return SourceBatch(new_batch)
