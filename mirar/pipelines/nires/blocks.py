"""
Script containing the various
:class:`~mirar.processors.base_processor.BaseProcessor`
lists which are used to build configurations for the
:class:`~mirar.pipelines.sedmv2.sedmv2_pipeline.SEDMv2Pipeline`.
"""

from mirar.paths import BASE_NAME_KEY, OBSCLASS_KEY, core_fields
from mirar.pipelines.nires.config import sextractor_astrometry_config
from mirar.pipelines.nires.load_nires_image import load_raw_nires_image
from mirar.processors.astromatic import Sextractor
from mirar.processors.astromatic.sextractor.background_subtractor import (
    SextractorBkgSubtractor,
)
from mirar.processors.csvlog import CSVLog
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.utils import (
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageSaver,
    ImageSelector,
)

load_raw = [
    ImageLoader(input_sub_dir="raw", load_image=load_raw_nires_image),
    ImageDebatcher(),
]

build_log = [  # pylint: disable=duplicate-code
    ImageBatcher(
        [
            "MJD-OBS",
        ]
    ),
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

flat_calibrate = [
    ImageSelector(("EXPTIME", "30.0")),
    ImageDebatcher(),
    ImageBatcher(
        [
            "FILTER",
        ]
    ),
    SkyFlatCalibrator(
        cache_sub_dir="calibration_flats",
        cache_image_name_header_keys=["FILTER", "TARGET"],
    ),
    ImageSaver(output_dir_name="skyflatcal"),
    ImageDebatcher(),
    ImageBatcher(
        [
            "MJD-OBS",
        ]
    ),
    Sextractor(
        **sextractor_astrometry_config,
        write_regions_bool=True,
        output_sub_dir="skysub",
        checkimage_type=["-BACKGROUND"],
    ),
    SextractorBkgSubtractor(),
    ImageSaver(output_dir_name="skysub"),
]
