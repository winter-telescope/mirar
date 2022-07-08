import logging
import astropy.io.fits
import numpy as np
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import base_name_key


logger = logging.getLogger(__name__)


class HeaderReader(BaseImageProcessor):

    base_key = "header_reader"

    def __init__(
            self,
            *keys,
    ):
        super().__init__()
        self.keys = keys

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for header in headers:

            msg = f"{header[base_name_key]} \n"
            for key in self.keys:
                msg += f"{key}: {header[key]} \n"
            logger.info(msg)

        return images, headers
