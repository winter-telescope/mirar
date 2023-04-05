"""
Script containing the various
:class:`~winterdrp.processors.base_processor.BaseProcessor`
lists which are used to build configurations for the
:class:`~winterdrp.pipelines.sedmv2.sedmv2_pipeline.SEDMv2Pipeline`.
"""
from winterdrp.downloader.get_test_data import get_test_data_dir
from winterdrp.paths import BASE_NAME_KEY, core_fields
from winterdrp.pipelines.sedmv2.config import (
    psfex_config_path,
    sedmv2_cal_requirements,
    sextractor_astrometry_config,
    sextractor_photometry_config,
    swarp_config_path,
)
from winterdrp.pipelines.sedmv2.generator import (
    sedmv2_photometric_catalog_generator,
    sedmv2_reference_image_generator,
    sedmv2_reference_image_resampler,
    sedmv2_reference_psfex,
    sedmv2_reference_sextractor,
)
from winterdrp.pipelines.sedmv2.load_sedmv2_image import load_raw_sedmv2_image
from winterdrp.processors import BiasCalibrator, FlatCalibrator
from winterdrp.processors.anet import AstrometryNet
from winterdrp.processors.astromatic import PSFex, Sextractor, Swarp
from winterdrp.processors.csvlog import CSVLog
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.reference import Reference
from winterdrp.processors.utils import (
    ImageBatcher,
    ImageLoader,
    ImageSaver,
    ImageSelector,
)
from winterdrp.processors.utils.cal_hunter import CalHunter
from winterdrp.processors.utils.header_annotate import HeaderEditor
from winterdrp.processors.utils.simulate_realtime import RealtimeImageSimulator
from winterdrp.processors.zogy.zogy import (
    ZOGY,
    ZOGYPrepare,
    default_sedmv2_catalog_purifier,
)

load_raw = [
    ImageLoader(load_image=load_raw_sedmv2_image),
]

load_test = [
    ImageLoader(
        input_img_dir=get_test_data_dir(),
        input_sub_dir="raw",
        load_image=load_raw_sedmv2_image,
    ),
]

sim_realtime = [
    RealtimeImageSimulator(
        input_img_dir=get_test_data_dir(),
        input_img_names=[
            "sedmv2/20220402/raw/SEDMV2_20220402_193104_Camera0.fits",
            "sedmv2/20220402/raw/SEDMV2_20220402_214324_Camera0.fits",
        ],
        output_dir=get_test_data_dir(),
        output_dir_name="raw",
    )
]

cal_hunter = [
    CalHunter(load_image=load_raw_sedmv2_image, requirements=sedmv2_cal_requirements),
]

build_log = [  # pylint: disable=duplicate-code
    CSVLog(
        export_keys=[
            "UTC",
            "FIELDID",
            "FILTERID",
            "EXPTIME",
            "OBSTYPE",
            "RA",
            "DEC",
            "TARGTYPE",
            "PROGID",
            "PROGPI",
            BASE_NAME_KEY,
        ]
        + core_fields
    ),
]  # pylint: disable=duplicate-code

process_raw = [
    BiasCalibrator(),
    ImageSelector(("OBSTYPE", ["FLAT", "SCIENCE"])),
    ImageBatcher(split_key="filter"),
    FlatCalibrator(),
    ImageBatcher(split_key=BASE_NAME_KEY),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),  # pylint: disable=duplicate-code
    ImageSaver(output_dir_name="detrend", write_mask=True),
    # AutoAstrometry(pa=0, inv=True, pixel_scale=SEDMV2_PIXEL_SCALE),
    # AstrometryNet(scale_upper=0.1667, scale_lower=0.0333, scale_units='degw'),
    AstrometryNet(scale_bounds=(0.1667, 0.0333), scale_units="degw", downsample=2),
    ImageSaver(output_dir_name="a-net"),
    # ImageSaver(
    #    output_dir_name="detrend", write_mask=True
    # ),  # pylint: disable=duplicate-code
    Sextractor(
        output_sub_dir="sextractor",
        checkimage_name=None,
        checkimage_type=None,
        **sextractor_astrometry_config
    ),
    Swarp(
        swarp_config_path=swarp_config_path,
        include_scamp=False,
    ),
    ImageSaver(
        output_dir_name="processed", write_mask=True
    ),  # pylint: disable=duplicate-code
    Sextractor(
        output_sub_dir="photprocess",
        checkimage_type="BACKGROUND_RMS",
        **sextractor_photometry_config
    ),  # pylint: disable=duplicate-code
    PhotCalibrator(ref_catalog_generator=sedmv2_photometric_catalog_generator),
    ImageSaver(
        output_dir_name="processed",
        write_mask=True,
    ),
    HeaderEditor(edit_keys="procflag", values=1),
]

subtract = [
    ImageBatcher(split_key=BASE_NAME_KEY),
    ImageSelector(("OBSTYPE", "SCIENCE")),
    Reference(
        ref_image_generator=sedmv2_reference_image_generator,
        ref_psfex=sedmv2_reference_psfex,
        sextractor=sedmv2_reference_sextractor,
        swarp_resampler=sedmv2_reference_image_resampler,  # pylint: disable=duplicate-code
    ),
    Sextractor(
        output_sub_dir="subtract",
        cache=False,
        write_regions_bool=True,
        **sextractor_photometry_config
    ),
    PSFex(config_path=psfex_config_path, output_sub_dir="subtract", norm_fits=True),
    ImageSaver(output_dir_name="ref"),
    ZOGYPrepare(
        output_sub_dir="subtract",
        sci_zp_header_key="ZP_AUTO",
        catalog_purifier=default_sedmv2_catalog_purifier,
    ),
    ZOGY(output_sub_dir="subtract"),
]

# imsub = subtract + export_diff_to_db + extract_candidates
imsub = subtract
