import logging
from winterdrp.processors.base_processor import BaseProcessor
import numpy as np
from astropy.io import fits
from collections import Callable
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor

logger = logging.getLogger(__name__)


class FilterCandidates(BaseProcessor):

    def __init__(self,
                 *args,
                 **kwargs):
        super(FilterCandidates, self).__init__(*args, **kwargs)
        pass

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[fits.Header],
    ) -> tuple[list[np.ndarray], list[fits.Header]]:
        pass


class DetectCandidates(BaseProcessor):

    def __init__(self,
                 scorr_image_sextracor: Callable[Sextractor],
                 *args,
                 **kwargs):
        super(DetectCandidates, self).__init__(*args, **kwargs)
        pass

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[fits.Header],
    ) -> tuple[list[np.ndarray], list[fits.Header]]:
        pass
