import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.processors.autoastrometry.autoastrometry import run_autoastrometry_single
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)


class AutoAstrometry(BaseProcessor):

    base_key = "autoastrometry"

    def __init__(
            self,
            temp_output_sub_dir: str = "autoastrometry",
            write_crosscheck_files: bool = False,
            *args,
            **kwargs
    ):
        super(AutoAstrometry, self).__init__(*args, **kwargs)

        self.temp_output_sub_dir = temp_output_sub_dir
        self.write_crosscheck_files = write_crosscheck_files

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        sextractor_out_dir = get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

        try:
            os.makedirs(sextractor_out_dir)
        except OSError:
            pass

        for i, data in enumerate(images):
            header = headers[i]

            temp_path = os.path.join(sextractor_out_dir, header["BASENAME"])

            self.save_fits(data, header, temp_path)

            run_autoastrometry_single(
                img_path=temp_path,
                output_dir=sextractor_out_dir,
                write_crosscheck_files=self.write_crosscheck_files,
                overwrite=True
            )

            # Load up temp path image.header, then delete
            img, header = self.open_fits(temp_path)
            images[i] = img
            headers[i] = header
            os.remove(temp_path)
            logger.info(f"Loaded updated header, and deleted temporary file {temp_path}")

        return images, headers
