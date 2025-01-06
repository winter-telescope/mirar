"""
Script containing the various
:class:`~mirar.processors.base_processor.BaseProcessor`
lists which are used to build configurations for the
:class:`~mirar.pipelines.gmos.gmos_pipeline.GMOSPipeline`.

Note: this pipeline works on GMOS data that has already been processed by DRAGONS
"""

# pylint: disable=duplicate-code

from mirar.paths import BASE_NAME_KEY, OBSCLASS_KEY, TARGET_KEY, core_fields
from mirar.pipelines.gmos.config import (
    psfex_sci_config_path,
    scamp_path,
    sextractor_astrometry_config,
    sextractor_photometry_config,
    swarp_config_path,
)
from mirar.pipelines.gmos.config.constants import GMOS_PIXEL_SCALE
from mirar.pipelines.gmos.generator import (
    gmos_astrometric_catalog_generator,
    gmos_photometric_catalog_generator,
    gmos_reference_image_generator,
    gmos_reference_image_resampler,
    gmos_reference_psfex,
    gmos_reference_sextractor,
    gmos_zogy_catalogs_purifier,
)
from mirar.pipelines.gmos.gmos_mask import generate_gmos_mask
from mirar.pipelines.gmos.load_gmos_image import load_detrended_gmos_image
from mirar.processors.astromatic import PSFex, Scamp, Sextractor
from mirar.processors.astromatic.swarp import Swarp
from mirar.processors.astrometry.anet import AstrometryNet
from mirar.processors.astrometry.autoastrometry import AutoAstrometry
from mirar.processors.csvlog import CSVLog
from mirar.processors.mask import MaskPixelsFromFunction
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.processors.photometry import AperturePhotometry, PSFPhotometry
from mirar.processors.reference import ProcessReference
from mirar.processors.sources import (
    CSVExporter,
    ForcedPhotometryDetector,
    ImageUpdater,
    ParquetWriter,
)
from mirar.processors.sources.utils import RegionsWriter
from mirar.processors.utils import ImageBatcher, ImageLoader, ImageRebatcher, ImageSaver
from mirar.processors.zogy.zogy import ZOGY, ZOGYPrepare

load_raw = [
    ImageLoader(input_sub_dir="detrend", load_image=load_detrended_gmos_image),
    ImageBatcher(BASE_NAME_KEY),
    MaskPixelsFromFunction(mask_function=generate_gmos_mask),
    ImageSaver(output_dir_name="mask"),
]

build_log = [  # pylint: disable=duplicate-code
    CSVLog(
        export_keys=[
            TARGET_KEY,
            "OBJRA",
            "OBJDEC",
            "DATE-OBS",
            "FILTER",
            OBSCLASS_KEY,
            BASE_NAME_KEY,
        ]
        + core_fields
    ),
]  # pylint: disable=duplicate-code
calibrate = [
    ImageBatcher(split_key=BASE_NAME_KEY),
    AutoAstrometry(),
    Sextractor(
        output_sub_dir="sextractor",
        **sextractor_astrometry_config,
    ),
    Scamp(
        ref_catalog_generator=gmos_astrometric_catalog_generator,
        scamp_config_path=scamp_path,
        cache=False,
    ),
    ImageRebatcher(split_key=["target", "filter"]),
    ImageSaver(output_dir_name="astrometry"),
    Swarp(
        swarp_config_path=swarp_config_path,
        include_scamp=True,
    ),
    AstrometryNet(
        output_sub_dir="anet",
        timeout=120,
        use_sextractor=True,
    ),
    Sextractor(
        output_sub_dir="photprocess",
        checkimage_type="BACKGROUND_RMS",
        **sextractor_photometry_config,
    ),
    PhotCalibrator(
        ref_catalog_generator=gmos_photometric_catalog_generator,
    ),
    ImageSaver(
        output_dir_name="processed",
    ),
]

reduce = build_log + calibrate

subtract = [
    ProcessReference(
        ref_image_generator=gmos_reference_image_generator,
        ref_psfex=gmos_reference_psfex,
        sextractor=gmos_reference_sextractor,
        swarp_resampler=gmos_reference_image_resampler,  # pylint: disable=duplicate-code
        temp_output_subtract_dir="subtract",
    ),
    Sextractor(
        output_sub_dir="subtract",
        cache=False,
        write_regions_bool=False,
        **sextractor_photometry_config,
    ),
    PSFex(
        config_path=psfex_sci_config_path,
        output_sub_dir="subtract",
        # norm_fits=True,
    ),
    ZOGYPrepare(
        output_sub_dir="zogy",
        catalog_purifier=gmos_zogy_catalogs_purifier,
    ),
    ZOGY(output_sub_dir="zogy"),
    ImageSaver(output_dir_name="diff"),
    ForcedPhotometryDetector(ra_header_key="OBJRA", dec_header_key="OBJDEC"),
    RegionsWriter(output_dir_name="diff"),
    AperturePhotometry(
        aper_diameters=[
            2 / GMOS_PIXEL_SCALE,
            3 / GMOS_PIXEL_SCALE,
            4 / GMOS_PIXEL_SCALE,
            5 / GMOS_PIXEL_SCALE,
            10 / GMOS_PIXEL_SCALE,
        ],
        bkg_in_diameters=[
            2.5 / GMOS_PIXEL_SCALE,
            3.5 / GMOS_PIXEL_SCALE,
            4.5 / GMOS_PIXEL_SCALE,
            5.5 / GMOS_PIXEL_SCALE,
            10.5 / GMOS_PIXEL_SCALE,
        ],
        bkg_out_diameters=[
            5.5 / GMOS_PIXEL_SCALE,
            8.6 / GMOS_PIXEL_SCALE,
            9.5 / GMOS_PIXEL_SCALE,
            10.6 / GMOS_PIXEL_SCALE,
            15.6 / GMOS_PIXEL_SCALE,
        ],
        col_suffix_list=["2", "3", "4", "5", "10"],
        phot_cutout_half_size=100,
        zp_key="ZP_AUTO",
    ),
    PSFPhotometry(),
    ParquetWriter(output_dir_name="sources"),
    CSVExporter(output_dir_name="sources"),
    ImageUpdater(modify_dir_name="diff"),
]
