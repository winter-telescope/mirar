from mirar.paths import BASE_NAME_KEY, EXPTIME_KEY, OBSCLASS_KEY, TARGET_KEY
from mirar.pipelines.spring.generator import (
    spring_anet_sextractor_config_path_generator,
)
from mirar.pipelines.spring.load_spring_image import load_raw_spring_image
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.flat import SkyFlatCalibrator
from mirar.processors.utils import (
    ImageLoader,
    ImageRebatcher,
    ImageSaver,
    ImageSelector,
)

load_raw = [
    ImageLoader(input_sub_dir="raw", load_image=load_raw_spring_image),
]

csvlog = [
    ImageRebatcher(BASE_NAME_KEY),
    CSVLog(
        export_keys=[
            # Time / Identity
            "UTCISO",
            "DATE-OBS",
            "MJD-OBS",
            BASE_NAME_KEY,
            "FILENAME",
            # Targeting
            TARGET_KEY,
            "TARGNAME",
            "FIELDID",
            # Pointing
            "RADEG",
            "DECDEG",
            "AZIMUTH",
            "ALTITUDE",
            "AIRMASS",
            # Instrument / exposure
            "INSTRUME",
            "TELESCOP",
            "OBSERVAT",
            "FILTER",
            EXPTIME_KEY,
            OBSCLASS_KEY,
            # Other keys
            "FOCPOS",
            "TMP_CUR",
            "TMP_SET",
            "READOUTM",
            "MEDCOUNT",
        ]
    ),
    ImageRebatcher(BASE_NAME_KEY),
]

dark_calibrate = [
    ImageRebatcher([EXPTIME_KEY]),
    DarkCalibrator(
        cache_sub_dir="calibration_darks",
        cache_image_name_header_keys=[EXPTIME_KEY],
    ),
    ImageRebatcher(BASE_NAME_KEY),
    ImageSaver(output_dir_name="darkcal"),
    ImageSelector((OBSCLASS_KEY, ["science"])),
]

flat_calibrate = [
    ImageRebatcher([TARGET_KEY]),
    SkyFlatCalibrator(
        cache_sub_dir="sky_dither_flats",
        cache_image_name_header_keys=["FILTER", TARGET_KEY],
        flat_mode="median",
        try_load_cache=False,
    ),
    ImageSaver(output_dir_name="skyflatcal"),
]

astrometry = [
    ImageRebatcher(BASE_NAME_KEY),
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[0.3, 0.5],
        scale_units="app",
        use_sextractor=True,
        # parity="neg",
        search_radius_deg=0.15,
        sextractor_config_path=spring_anet_sextractor_config_path_generator,
        use_weight=True,
        timeout=120,
        cache=False,
        no_tweak=True,
    ),
    ImageSaver("post_astrometry"),
]
