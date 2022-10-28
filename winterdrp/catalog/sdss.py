import logging
from pathlib import Path

import numpy as np
from pydl.pydlutils import mangle
import requests

from winterdrp.catalog.base_catalog import VizierCatalog
from winterdrp.errors import ProcessorError
from winterdrp.paths import base_output_dir


logger = logging.getLogger(__name__)


class NotInSDSSError(ProcessorError):
    pass


SDSS_RELEASE = "dr16"

SDSS_COVERAGE_PATH = Path(base_output_dir).joinpath(f"SDSS/{SDSS_RELEASE}_window.ply")

SDSS_COVERAGE_URL = f"https://data.sdss.org/sas/{SDSS_RELEASE}/sdss/tiling/final/internal/tmp_window_unif.ply"


def get_sdss_coverage(
) -> mangle.PolygonList:

    if not SDSS_COVERAGE_PATH.parent.exists():
        SDSS_COVERAGE_PATH.parent.mkdir(parents=True)

    if not SDSS_COVERAGE_PATH.is_file():
        logger.info(f"No coverage found. Downloading SDSS coverage map from {SDSS_COVERAGE_URL}")
        print(SDSS_COVERAGE_PATH)
        with open(str(SDSS_COVERAGE_PATH), "wb+") as f:
            r = requests.get(SDSS_COVERAGE_URL, allow_redirects=True)
            f.write(r.content)

    sdss_coverage = mangle.read_mangle_polygons(str(SDSS_COVERAGE_PATH))
    return sdss_coverage


def in_sdss(
        ra_deg: float,
        dec_deg: float
) -> bool:
    sdss_coverage = get_sdss_coverage()
    return mangle.is_in_window(sdss_coverage, np.array([[ra_deg, dec_deg]]))[0]


class SDSS(VizierCatalog):

    catalog_vizier_code = "V/154"
    abbreviation = "sdss"

    ra_key = "RA_ICRS"
    dec_key = "DE_ICRS"

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        if not in_sdss(ra_deg, dec_deg):
            err = f"Querying for SDSS sources, but the field ({ra_deg}, {dec_deg}) was not observed in SDSS."
            logger.error(err)
            raise NotInSDSSError(err)
