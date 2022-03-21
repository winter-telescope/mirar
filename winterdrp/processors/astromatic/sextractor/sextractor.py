import os
import numpy as np
import logging
import astropy.io.fits
from winterdrp.paths import astrometry_output_dir
from winterdrp.processors.astromatic.sextractor.sourceextractor import run_sextractor_single, default_config
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)


class Sextractor(BaseProcessor):

    base_key = "sextractor"

    def __init__(
            self,
            output_sub_dir: str,
            config: str = default_config,
            *args,
            **kwargs
    ):
        super(Sextractor, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir
        self.config = config

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        sextractor_out_dir = get_output_dir(self.output_sub_dir, self.night_sub_dir)

        try:
            os.makedirs(sextractor_out_dir)
        except OSError:
            pass

        for i, data in enumerate(images):
            header = headers[i]

            temp_path = os.path.join(sextractor_out_dir, header["BASENAME"])

            self.save_fits(data, header, temp_path)

            run_sextractor_single(
                img=temp_path,
                config=self.config,
                output_dir=sextractor_out_dir,
            )

            # Load up temp path image.header, then delete
            img, header = self.open_fits(temp_path)
            images[i] = img
            headers[i] = header
            os.remove(temp_path)
            logger.info(f"Loaded updated header, and deleted temporary file {temp_path}")

        return images, headers

    # def _apply_to_images(
    #         self,
    #         images: list,
    #         headers: list,
    #         sub_dir: str = ""
    # ) -> (list, list):
    #
    #     # Try making output directory, unless it exists
    #
    #     output_dir = astrometry_output_dir(sub_dir)
    #
    #     try:
    #         os.makedirs(output_dir)
    #     except OSError:
    #         pass
    #
    #     for header in list(headers):
    #
    #         # First run Sextractor
    #
    #         run_sextractor(
    #             header[latest_save_key],
    #             config=self.config,
    #             output_dir=output_dir,
    #         )
    #
    #     return images, headers
