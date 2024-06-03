"""
Module containing reference-related functions for the winter pipeline.
"""

import logging
import os

from mirar.data import Image, ImageBatch
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions import select_from_table
from mirar.paths import REF_CAT_PATH_KEY, get_output_dir
from mirar.pipelines.winter.config import (
    ref_psfex_path,
    sextractor_astrometry_config,
    sextractor_reference_config,
    sextractor_reference_psf_phot_config,
    swarp_config_path,
)
from mirar.pipelines.winter.constants import winter_filters_map
from mirar.pipelines.winter.generator.utils import winter_ref_catalog_namer
from mirar.pipelines.winter.models import (
    DEFAULT_FIELD,
    RefComponent,
    RefQuery,
    RefStack,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.processors.split import SUB_ID_KEY
from mirar.references import PS1Ref
from mirar.references.local import RefFromPath
from mirar.references.wfcam.wfcam_query import WFAUQuery
from mirar.references.wfcam.wfcam_stack import WFCAMStackedRef

logger = logging.getLogger(__name__)


def winter_astrometric_ref_catalog_namer(batch: ImageBatch) -> ImageBatch:
    """
    Function to name the reference catalog to use for WINTER astrometry,
    and add this to header

    :param batch: ImageBatch
    :return: ImageBatch
    """
    winter_reference_catalog_dir = get_output_dir(
        dir_root="astrometric", sub_dir="winter/reference_catalogs"
    )
    for ind, image in enumerate(batch):
        image[REF_CAT_PATH_KEY] = winter_ref_catalog_namer(
            image, winter_reference_catalog_dir
        ).as_posix()
        batch[ind] = image
    return batch


# Swarp generators
def winter_reference_image_resampler_for_zogy(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(
        swarp_config_path=swarp_config_path, cache=True, subtract_bkg=False, **kwargs
    )


def winter_wfau_component_image_stacker(**kwargs) -> Swarp:
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


def winter_reference_sextractor(
    output_sub_dir: str,
) -> Sextractor:
    """
    Returns a Sextractor processor for WINTER reference images

    :param output_sub_dir: str
    :return: Sextractor processor
    """
    return Sextractor(
        **sextractor_reference_config,
        output_sub_dir=output_sub_dir,
        cache=True,
    )


def winter_reference_psf_phot_sextractor(output_sub_dir: str) -> Sextractor:
    """
    Returns a Sextractor processor for WINTER reference images

    :param output_sub_dir: str
    :return: Sextractor processor
    """
    return Sextractor(
        **sextractor_reference_psf_phot_config,
        output_sub_dir=output_sub_dir,
        cache=False,
        use_psfex=True,
    )


def winter_reference_psfex(output_sub_dir: str, norm_fits: bool) -> PSFex:
    """
    Returns a PSFEx processor for WINTER

    :param output_sub_dir: str
    :param norm_fits: bool
    """
    return PSFex(
        config_path=ref_psfex_path,
        output_sub_dir=output_sub_dir,
        norm_fits=norm_fits,
    )


def ref_sextractor(image: Image):
    """
    Generates a sextractor instance for reference images to get photometry

    :param image: Image
    :return: Sextractor
    """
    logger.debug(image)
    return Sextractor(
        output_sub_dir="phot",
        **sextractor_astrometry_config,
        write_regions_bool=True,
        cache=False,
    )


def winter_reference_stack_annotator(stacked_image: Image, image: Image) -> Image:
    """
    Generates a stack id for WINTER reference images

    :param stacked_image: Image
    :param image: Image
    """
    stackid = (
        f"{str(image.header['FIELDID']).rjust(5, '0')}"
        f"{str(image.header[SUB_ID_KEY]).rjust(2, '0')}"
        f"{str(winter_filters_map[image.header['FILTER']])}"
    )
    stacked_image["STACKID"] = int(stackid)
    stacked_image["FIELDID"] = image.header["FIELDID"]
    stacked_image[SUB_ID_KEY] = image.header[SUB_ID_KEY]
    return stacked_image


def winter_reference_generator(image: Image):
    """
    Gets a reference image generator for the winter data

    :param image: Image
    :return: Reference image generator
    """
    components_image_dir = get_output_dir(
        dir_root="components", sub_dir="winter/references"
    )
    components_image_dir.mkdir(parents=True, exist_ok=True)

    filtername = image["FILTER"]

    if filtername not in ["Y", "J", "H"]:
        raise ValueError(f"Filter {filtername} not recognized for WINTER")

    # TODO if in_ukirt and in_vista, different processing
    fieldid = int(image["FIELDID"])
    subdetid = int(image[SUB_ID_KEY])
    logger.debug(f"Fieldid: {fieldid}, subdetid: {subdetid}")

    cache_ref_stack = False
    if filtername in ["J", "H"]:
        if fieldid != DEFAULT_FIELD:
            cache_ref_stack = True
            constraints = DBQueryConstraints(
                columns=["fieldid", SUB_ID_KEY.lower()],
                accepted_values=[fieldid, subdetid],
            )

            db_results = select_from_table(
                db_constraints=constraints,
                sql_table=RefStack.sql_model,
                output_columns=["savepath"],
            )

            if len(db_results) > 0:
                savepath = db_results["savepath"].iloc[0]
                if os.path.exists(savepath):
                    logger.debug(f"Found reference image in database: {savepath}")
                    return RefFromPath(path=savepath, filter_name=filtername)

        skip_online_query = filtername == "H"

        wfcam_query = WFAUQuery(
            num_query_points=16,
            filter_name=filtername,
            use_db_for_component_queries=True,
            components_db_table=RefComponent,
            query_db_table=RefQuery,
            skip_online_query=skip_online_query,
            component_image_subdir="winter/references/components",
        )

        return WFCAMStackedRef(
            filter_name=filtername,
            wfcam_query=wfcam_query,
            image_resampler_generator=winter_wfau_component_image_stacker,
            write_stacked_image=cache_ref_stack,
            write_stack_sub_dir="winter/references/ref_stacks",
            write_stack_to_db=cache_ref_stack,
            stacks_db_table=RefStack,
            component_image_sub_dir="components",
            references_base_subdir_name="winter/references",
            stack_image_annotator=winter_reference_stack_annotator,
        )

    assert filtername == "Y", f"Filter {filtername} not recognized for WINTER"

    # Use PS1 references for Y-band
    logger.debug("Will query reference image from PS1")
    return PS1Ref(filter_name=filtername)


def winter_photometric_ref_catalog_namer(batch: ImageBatch) -> ImageBatch:
    """
    Function to name the reference catalog to use for WINTER astrometry
    """
    winter_reference_catalog_dir = get_output_dir(
        dir_root="photometric", sub_dir="winter/reference_catalogs"
    )
    winter_reference_catalog_dir.mkdir(exist_ok=True, parents=True)
    for ind, image in enumerate(batch):
        image[REF_CAT_PATH_KEY] = winter_ref_catalog_namer(
            image, winter_reference_catalog_dir
        ).as_posix()
        batch[ind] = image
    return batch
