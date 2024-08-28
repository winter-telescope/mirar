"""
Script containing the various
:class:`~mirar.processors.base_processor.BaseProcessor`
lists which are used to build configurations for the
:class:`~mirar.pipelines.neowise.neowise_pipeline.NEOWISEPipeline`.
"""

from mirar.paths import (
    BASE_NAME_KEY,
    GAIN_KEY,
    OBSCLASS_KEY,
    ZP_KEY,
    ZP_STD_KEY,
    core_fields,
)
from mirar.pipelines.neowise.config import (
    NEOWISE_PIXEL_SCALE,
    psfex_sci_config_path,
    sextractor_photometry_config,
)

# from mirar.pipelines.neowise.config.constants import GIT_PIXEL_SCALE
from mirar.pipelines.neowise.generator import (
    neowise_reference_image_generator,
    neowise_reference_image_resampler,
    neowise_reference_psfex,
    neowise_reference_sextractor,
    neowise_zogy_catalogs_purifier,
)
from mirar.pipelines.neowise.load_neowise_image import load_stacked_neowise_image
from mirar.processors.astromatic import PSFex, Sextractor
from mirar.processors.csvlog import CSVLog
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import ProcessReference
from mirar.processors.sources import ForcedPhotometryDetector, SourceWriter
from mirar.processors.utils import ImageBatcher, ImageLoader, ImageSaver, ImageSelector
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_stack_neowise = [
    ImageLoader(input_sub_dir="stack", load_image=load_stacked_neowise_image),
    ImageBatcher(split_key=BASE_NAME_KEY),
]

load_ref_neowise = [
    ImageLoader(input_sub_dir="icore_reference", load_image=load_stacked_neowise_image)
]
save_image = [ImageSaver(output_dir_name="references")]
subtract = [
    ImageBatcher(split_key=BASE_NAME_KEY),
    ImageSelector((OBSCLASS_KEY, "science")),
    ProcessReference(
        ref_image_generator=neowise_reference_image_generator,
        ref_psfex=neowise_reference_psfex,
        sextractor=neowise_reference_sextractor,
        swarp_resampler=neowise_reference_image_resampler,
        # pylint: disable=duplicate-code
        temp_output_subtract_dir="subtract",
    ),
    Sextractor(
        output_sub_dir="subtract",
        cache=False,
        write_regions_bool=True,
        **sextractor_photometry_config,
    ),
    PSFex(
        config_path=psfex_sci_config_path,
        output_sub_dir="subtract",
        norm_fits=True,
    ),
    ImageSaver(output_dir_name="ref"),
    ZOGYPrepare(
        output_sub_dir="subtract",
        catalog_purifier=neowise_zogy_catalogs_purifier,
    ),
    ZOGY(output_sub_dir="subtract"),
    ImageSaver(output_dir_name="diff"),
]

build_log = [  # pylint: disable=duplicate-code
    CSVLog(
        export_keys=[
            "DATE-OBS",
            "FILTER",
            OBSCLASS_KEY,
            BASE_NAME_KEY,
        ]
        + core_fields
    ),
]  # pylint: disable=duplicate-code

image_photometry = [  # imported from wirc/blocks.py
    ImageSaver(output_dir_name="for_photometry"),
    ForcedPhotometryDetector(ra_header_key="OBJRAD", dec_header_key="OBJDECD"),
    AperturePhotometry(
        aper_diameters=[
            2 / NEOWISE_PIXEL_SCALE,
            3 / NEOWISE_PIXEL_SCALE,
            4 / NEOWISE_PIXEL_SCALE,
            5 / NEOWISE_PIXEL_SCALE,
            10 / NEOWISE_PIXEL_SCALE,
        ],
        bkg_in_diameters=[
            2.5 / NEOWISE_PIXEL_SCALE,
            3.5 / NEOWISE_PIXEL_SCALE,
            4.5 / NEOWISE_PIXEL_SCALE,
            5.5 / NEOWISE_PIXEL_SCALE,
            10.5 / NEOWISE_PIXEL_SCALE,
        ],
        bkg_out_diameters=[
            5.5 / NEOWISE_PIXEL_SCALE,
            8.6 / NEOWISE_PIXEL_SCALE,
            9.5 / NEOWISE_PIXEL_SCALE,
            10.6 / NEOWISE_PIXEL_SCALE,
            15.6 / NEOWISE_PIXEL_SCALE,
        ],
        col_suffix_list=["2", "3", "4", "5", "10"],
        phot_cutout_half_size=100,
        zp_key="ZP",
    ),
    PSFPhotometry(),
    SourceWriter(output_dir_name="photometry_table"),
]

candidate_photometry = [  # imported from wirc/blocks.py
    AperturePhotometry(
        aper_diameters=[16, 70],
        phot_cutout_half_size=100,
        bkg_in_diameters=[25, 90],
        bkg_out_diameters=[40, 100],
        col_suffix_list=["", "big"],
    ),
    PSFPhotometry(),
]
