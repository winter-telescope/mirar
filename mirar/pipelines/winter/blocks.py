from mirar.pipelines.winter.generator import (
    winter_mask_path,
    winter_reference_generator,
)
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.mask import MaskPixelsFromPath
from mirar.processors.reference import GetReferenceImage
from mirar.processors.sky import NightSkyMedianCalibrator, SkyFlatCalibrator
from mirar.processors.split import SplitImage
from mirar.processors.utils import ImageDebatcher, ImageSaver
from mirar.processors.utils.multi_ext_parser import MultiExtParser

refbuild = [
    ImageDebatcher(),
    GetReferenceImage(
        ref_image_generator=winter_reference_generator,
    ),
    ImageSaver(output_dir_name="stacked_ref"),
]

commissioning = [
    CSVLog(
        export_keys=[
            "OBJECT",
            "FILTER",
            "UTSHUT",
            "EXPTIME",
            "COADDS",
            "OBSTYPE",
            "OBSCLASS",
        ]
    ),
    MultiExtParser(input_sub_dir="raw/mef/"),
    SplitImage(),
    MaskPixelsFromPath(mask_path=winter_mask_path),
    DarkCalibrator(),
    SkyFlatCalibrator(),
    NightSkyMedianCalibrator(),
]
