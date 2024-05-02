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
from mirar.data import Image, SourceBatch
from mirar.pipelines.wirc.wirc_files import (
    psfex_path,
    sextractor_reference_config,
    wirc_file_dir,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.references.wirc import WIRCRef

logger = logging.getLogger(__name__)


def wirc_source_table_filter_annotator(source_table: SourceBatch) -> SourceBatch:
    """
    Function to remove bad candidates with None in sigmapsf, magpsf, magap, sigmagap,
    and update the source table with the keys required for the WIRC database
    :param source_table: source table
    :return: updated source table
    """

    new_batch = SourceBatch([])

    for source in source_table:
        src_df = source.get_data()

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

        source.set_data(src_df)
        new_batch.append(source)

    return new_batch


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
    if filter_name == "ks":
        return Gaia2Mass(
            min_mag=10, max_mag=20, search_radius_arcmin=10, filter_name="k"
        )
    return Gaia2Mass(
        min_mag=10, max_mag=20, search_radius_arcmin=10, filter_name=filter_name
    )


# def wirc_reference_image_generator(
#     image: Image,
#     images_directory: str = os.getenv("REF_IMG_DIR"),
# ) -> WIRCRef:
#     """
#     Function to match a new wirc image to a reference image directory
#
#     :param image: image
#     :param images_directory: ref image directory
#     :return: wirc ref
#     """
#     object_name = image["OBJECT"]
#     filter_name = image["FILTER"]
#     return WIRCRef(
#         object_name=object_name,
#         filter_name=filter_name,
#         images_directory_path=images_directory,
#     )

# def wirc_reference_image_generator(
#     image: Image,
# ):
#     def winter_reference_generator(image: Image):
#         """
#         Generates a reference image for the winter data
#         Args:
#             db_table: Database table to search for existing image
#             image: Image
#
#         Returns:
#
#         """
#         components_image_dir = get_output_dir(
#             dir_root="components", sub_dir="winter/references"
#         )
#         components_image_dir.mkdir(parents=True, exist_ok=True)
#
#         filtername = image["FILTER"]
#
#         if filtername not in ["Y", "J", "H"]:
#             raise ValueError(f"Filter {filtername} not recognized for WINTER")
#
#         # TODO if in_ukirt and in_vista, different processing
#         fieldid = int(image["FIELDID"])
#         subdetid = int(image[SUB_ID_KEY])
#         logger.debug(f"Fieldid: {fieldid}, subdetid: {subdetid}")
#
#         cache_ref_stack = False
#         if filtername in ["J", "H"]:
#             if fieldid != DEFAULT_FIELD:
#                 cache_ref_stack = True
#                 constraints = DBQueryConstraints(
#                     columns=["fieldid", SUB_ID_KEY.lower()],
#                     accepted_values=[fieldid, subdetid],
#                 )
#
#                 db_results = select_from_table(
#                     db_constraints=constraints,
#                     sql_table=RefStack.sql_model,
#                     output_columns=["savepath"],
#                 )
#
#                 if len(db_results) > 0:
#                     savepath = db_results["savepath"].iloc[0]
#                     if os.path.exists(savepath):
#                         logger.debug(f"Found reference image in database: {savepath}")
#                         return RefFromPath(path=savepath, filter_name=filtername)
#
#             skip_online_query = filtername == "H"
#
#             wfcam_query = WFAUQuery(
#                 num_query_points=16,
#                 filter_name=filtername,
#                 use_db_for_component_queries=True,
#                 components_db_table=RefComponent,
#                 query_db_table=RefQuery,
#                 skip_online_query=skip_online_query,
#                 component_image_subdir="winter/references/components",
#             )
#
#             return WFCAMStackedRef(
#                 filter_name=filtername,
#                 wfcam_query=wfcam_query,
#                 image_resampler_generator=winter_wfau_component_image_stacker,
#                 write_stacked_image=cache_ref_stack,
#                 write_stack_sub_dir="winter/references/ref_stacks",
#                 write_stack_to_db=cache_ref_stack,
#                 stacks_db_table=RefStack,
#                 component_image_sub_dir="components",
#                 references_base_subdir_name="winter/references",
#                 stack_image_annotator=winter_reference_stack_annotator,
#             )
#
#         assert filtername == "Y", f"Filter {filtername} not recognized for WINTER"
#
#         # Use PS1 references for Y-band
#         logger.debug("Will query reference image from PS1")
#         return PS1Ref(filter_name=filtername)


def wirc_reference_image_resampler(**kwargs) -> Swarp:
    """Returns a SWarp resampler for WIRC"""
    return Swarp(
        swarp_config_path=wirc_file_dir.joinpath("config.swarp"),
        cache=True,
        subtract_bkg=True,
        **kwargs
    )


def wirc_reference_sextractor(output_sub_dir: str) -> Sextractor:
    """Returns a Sextractor processor for WIRC reference images"""
    return Sextractor(
        **sextractor_reference_config, output_sub_dir=output_sub_dir, cache=True
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
