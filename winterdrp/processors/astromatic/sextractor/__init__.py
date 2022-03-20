import logging
import os

import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.processors.utils.image_saver import latest_save_key, ImageSaver
from winterdrp.processors.astromatic.sextractor.autoastrometry import run_autoastrometry
from winterdrp.paths import output_dir

logger = logging.getLogger(__name__)


class AutoAstrometry(BaseProcessor):

    base_key = "autoastrometry"

    def __init__(
            self,
            output_sub_dir: str,
            *args,
            **kwargs
    ):
        super(AutoAstrometry, self).__init__(*args, **kwargs)

        self.output_sub_dir = output_sub_dir

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        sextractor_out_dir = output_dir(self.output_sub_dir, self.night_sub_dir)

        try:
            os.makedirs(sextractor_out_dir)
        except OSError:
            pass

        for header in headers:
            run_autoastrometry(
                files=header[latest_save_key],
                output_dir=sextractor_out_dir
            )

            raise

        return images, headers

    def check_prerequisites(
            self,
    ):
        if not isinstance(self.preceding_steps[-1], ImageSaver):
            err = f"{self.__class__} requires {ImageSaver} to be run as the preceding step. " \
                  f"The following preceding steps were found: \n {self.preceding_steps}"
            logger.error(err)
            raise ValueError(err)
