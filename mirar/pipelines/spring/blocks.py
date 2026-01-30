from mirar.paths import BASE_NAME_KEY, EXPTIME_KEY, OBSCLASS_KEY, TARGET_KEY
from mirar.pipelines.spring.load_spring_image import load_raw_spring_image
from mirar.processors.csvlog import CSVLog
from mirar.processors.utils import ImageRebatcher
from mirar.processors.utils.image_loader import ImageLoader

load_raw = [ImageLoader(input_sub_dir="raw", load_image=load_raw_spring_image)]
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
        ]
    ),
    ImageRebatcher(BASE_NAME_KEY),
]
