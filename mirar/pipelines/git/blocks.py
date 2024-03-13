"""
Script containing the various
:class:`~mirar.processors.base_processor.BaseProcessor`
lists which are used to build configurations for the
:class:`~mirar.pipelines.sedmv2.sedmv2_pipeline.SEDMv2Pipeline`.
"""

from mirar.paths import (
    BASE_NAME_KEY,
    GAIN_KEY,
    OBSCLASS_KEY,
    ZP_KEY,
    ZP_STD_KEY,
    core_fields,
)
from mirar.pipelines.git.config import (
    psfex_sci_config_path,
    sextractor_astrometry_config,
    sextractor_photometry_config,
    swarp_config_path,
)
from mirar.pipelines.git.config.constants import GIT_PIXEL_SCALE
from mirar.pipelines.git.generator import (
    git_reference_image_generator,
    git_reference_image_resampler,
    git_reference_psfex,
    git_reference_sextractor,
    git_zogy_catalogs_purifier,
    lt_photometric_catalog_generator,
)
from mirar.pipelines.git.load_git_image import load_raw_git_image, load_raw_lt_image
from mirar.processors.astromatic import PSFex, Sextractor
from mirar.processors.astromatic.sextractor.background_subtractor import (
    SextractorBkgSubtractor,
)
from mirar.processors.astromatic.swarp import Swarp
from mirar.processors.cosmic_rays import LACosmicCleaner
from mirar.processors.csvlog import CSVLog
from mirar.processors.photcal import OutlierRejectionZPCalculator
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import ProcessReference
from mirar.processors.sources import ForcedPhotometryDetector
from mirar.processors.utils import (
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageSaver,
    ImageSelector,
)
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [
    ImageLoader(input_sub_dir="stack", load_image=load_raw_git_image),
    ImageSelector(("FILTER", "r")),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    ImageSaver(output_dir_name="skysub"),
    ImageDebatcher(),
    ImageBatcher(BASE_NAME_KEY),
    Swarp(
        center_ra=148.9996947,
        center_dec=69.6747746,
        center_type="MANUAL",
        x_imgpixsize=1500,
        y_imgpixsize=1500,
        swarp_config_path=swarp_config_path,
        include_scamp=False,
        propogate_headerlist=[ZP_KEY, ZP_STD_KEY],
    ),
    ImageSaver(output_dir_name="swarp_cut"),
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
    ImageSaver(output_dir_name="photometry"),
    ForcedPhotometryDetector(ra_header_key="OBJRAD", dec_header_key="OBJDECD"),
    AperturePhotometry(
        aper_diameters=[
            2 / GIT_PIXEL_SCALE,
            3 / GIT_PIXEL_SCALE,
            4 / GIT_PIXEL_SCALE,
            5 / GIT_PIXEL_SCALE,
            10 / GIT_PIXEL_SCALE,
        ],
        bkg_in_diameters=[
            2.5 / GIT_PIXEL_SCALE,
            3.5 / GIT_PIXEL_SCALE,
            4.5 / GIT_PIXEL_SCALE,
            5.5 / GIT_PIXEL_SCALE,
            10.5 / GIT_PIXEL_SCALE,
        ],
        bkg_out_diameters=[
            5.5 / GIT_PIXEL_SCALE,
            8.6 / GIT_PIXEL_SCALE,
            9.5 / GIT_PIXEL_SCALE,
            10.6 / GIT_PIXEL_SCALE,
            15.6 / GIT_PIXEL_SCALE,
        ],
        col_suffix_list=["2", "3", "4", "5", "10"],
        phot_cutout_half_size=100,
        zp_key="ZP_AUTO",
    ),
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


# split_stack = [
#     ImageDebatcher(),
#     ImageBatcher([BASE_NAME_KEY]),
#     SwarpImageSplitter(swarp_config_path=swarp_config_path, n_x=2, n_y=1),
#     ImageSaver(output_dir_name="split_stacks"),
# ]

# transients --

subtract = [
    ImageBatcher(split_key=BASE_NAME_KEY),
    ImageSelector((OBSCLASS_KEY, "science")),
    ProcessReference(
        ref_image_generator=git_reference_image_generator,
        ref_psfex=git_reference_psfex,
        sextractor=git_reference_sextractor,
        swarp_resampler=git_reference_image_resampler,  # pylint: disable=duplicate-code
        temp_output_subtract_dir="subtract_sdss",
    ),
    Sextractor(
        output_sub_dir="subtract_sdss",
        cache=False,
        write_regions_bool=True,
        **sextractor_photometry_config,
    ),
    PSFex(
        config_path=psfex_sci_config_path,
        output_sub_dir="subtract_sdss",
        norm_fits=True,
    ),
    ImageSaver(output_dir_name="ref"),
    ZOGYPrepare(
        output_sub_dir="subtract_sdss",
        catalog_purifier=git_zogy_catalogs_purifier,
    ),
    ZOGY(output_sub_dir="subtract_sdss"),
    ImageSaver(output_dir_name="diff_sdss"),
]

imsub = subtract  # + export_diff_to_db + extract_candidates


reduce_raw_lt = [
    ImageLoader(input_sub_dir="raw", load_image=load_raw_lt_image),
    ImageSelector(("FILTER", "r")),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    LACosmicCleaner(effective_gain_key=GAIN_KEY, readnoise=2),
    ImageSaver(output_dir_name="skysub"),
    ImageBatcher("FILTER"),
    Swarp(swarp_config_path=swarp_config_path, include_scamp=False),
    Sextractor(
        output_sub_dir="photometry",
        **sextractor_photometry_config,
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=lt_photometric_catalog_generator,
        zp_calculator=OutlierRejectionZPCalculator(
            outlier_rejection_threshold=[1.5, 2.0, 3.0]
        ),
    ),
    ImageSaver(output_dir_name="stack"),
]


load_stack_lt = [ImageLoader(input_sub_dir="stack")]
