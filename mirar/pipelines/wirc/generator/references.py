"""
Module for generating reference images for the WIRC pipeline
"""

import logging

from mirar.data import Image
from mirar.paths import get_output_dir
from mirar.pipelines.winter.generator.references import (
    winter_wfau_component_image_stacker,
)
from mirar.pipelines.wirc.wirc_files import (
    psfex_path,
    sextractor_reference_config,
    wirc_file_dir,
)
from mirar.processors.astromatic import PSFex, Sextractor, Swarp
from mirar.references.wfcam.wfcam_query import WFAUQuery
from mirar.references.wfcam.wfcam_stack import WFCAMStackedRef

logger = logging.getLogger(__name__)


def wirc_reference_generator(image: Image):
    """
    Gets a reference image generator for the wirc data

    :param image: Image
    :return: Reference image generator
    """
    components_image_dir = get_output_dir(
        dir_root="components", sub_dir="wirc/references"
    )
    components_image_dir.mkdir(parents=True, exist_ok=True)

    filtername = image["FILTER"]

    if filtername not in ["J", "H", "Ks"]:
        raise ValueError(f"Filter {filtername} not recognized for WINTER")

    cache_ref_stack = False

    wfcam_query = WFAUQuery(
        num_query_points=16,
        filter_name=filtername,
        use_db_for_component_queries=False,
        skip_online_query=False,
        component_image_subdir="wirc/references/components",
    )

    return WFCAMStackedRef(
        filter_name=filtername,
        wfcam_query=wfcam_query,
        image_resampler_generator=winter_wfau_component_image_stacker,
        write_stacked_image=cache_ref_stack,
        write_stack_sub_dir="wirc/references/ref_stacks",
        write_stack_to_db=cache_ref_stack,
        component_image_sub_dir="components",
        references_base_subdir_name="wirc/references",
        # stack_image_annotator=winter_reference_stack_annotator,
    )


def wirc_reference_image_resampler(**kwargs) -> Swarp:
    """Returns a SWarp resampler for WIRC"""
    return Swarp(
        swarp_config_path=wirc_file_dir.joinpath("config.swarp"),
        cache=True,
        subtract_bkg=True,
        **kwargs,
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
