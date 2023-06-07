import logging
from collections.abc import Callable

import astropy.table
from astropy.io import fits

from mirar.catalog.base_catalog import BaseCatalog
from mirar.processors.base_processor import BaseDataframeProcessor

logger = logging.getLogger(__name__)


class BrightStarCandidateMask(BaseDataframeProcessor):
    def __init__(
        self,
        ref_catalog_generator: Callable[[astropy.io.fits.Header], BaseCatalog],
        bright_star_thresh_mag: float = 14,
        *args,
        **kwargs
    ):
        super(BrightStarCandidateMask, self).__init__(*args, **kwargs)
        self.ref_catalog_generator = ref_catalog_generator
        self.bright_star_thresh_mag = bright_star_thresh_mag
