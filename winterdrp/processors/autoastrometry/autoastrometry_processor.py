import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.processors.autoastrometry.autoastrometry import run_autoastrometry_single
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)


class AutoAstrometry(BaseImageProcessor):

    base_key = "autoastrometry"

    def __init__(
            self,
            temp_output_sub_dir: str = "autoastrometry",
            write_crosscheck_files: bool = False,
            catalog: str = None,
            pixel_scale: float = None,
            inv: bool = False,
            pa: float = None,
            *args,
            **kwargs
    ):
        super(AutoAstrometry, self).__init__(*args, **kwargs)

        self.temp_output_sub_dir = temp_output_sub_dir
        self.write_crosscheck_files = write_crosscheck_files
        self.catalog = catalog
        self.pixel_scale = pixel_scale
        self.inv = inv
        self.pa = pa

    def __str__(self) -> str:
        return f"Processor to perform astrometric calibration."

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
                overwrite=True,
                catalog=self.catalog,
                pixel_scale=self.pixel_scale,
                inv=self.inv,
                pa=self.pa
            )

            # Load up temp path image.header, then delete
            img, header = self.open_fits(temp_path)

            images[i] = img
            headers[i] = header

            os.remove(temp_path)
            logger.info(f"Loaded updated header, and deleted temporary file {temp_path}")

        return images, headers
