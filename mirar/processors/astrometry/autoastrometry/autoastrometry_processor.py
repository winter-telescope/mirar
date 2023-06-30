"""
Module containing a processor to run automastrometry
"""
import logging
from typing import Optional

from mirar.data import ImageBatch
from mirar.paths import BASE_NAME_KEY, get_output_dir
from mirar.processors.astrometry.autoastrometry.autoastrometry import (
    run_autoastrometry_single,
)
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class AutoAstrometry(BaseImageProcessor):
    """Processor to run automastrometry"""

    base_key = "autoastrometry"

    def __init__(
        self,
        temp_output_sub_dir: str = "autoastrometry",
        write_crosscheck_files: bool = False,
        catalog: Optional[str] = None,
        pixel_scale: Optional[float] = None,
        inv: bool = False,
        pa: Optional[float] = None,
    ):
        super().__init__()

        self.temp_output_sub_dir = temp_output_sub_dir
        self.write_crosscheck_files = write_crosscheck_files
        self.catalog = catalog
        self.pixel_scale = pixel_scale
        self.inv = inv
        self.pa = pa

    def __str__(self) -> str:
        return "Processor to perform astrometric calibration."

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        sextractor_out_dir = get_output_dir(
            self.temp_output_sub_dir, self.night_sub_dir
        )

        sextractor_out_dir.mkdir(parents=True, exist_ok=True)

        for i, image in enumerate(batch):
            temp_path = sextractor_out_dir.joinpath(image[BASE_NAME_KEY])
            self.save_fits(image, temp_path)

            run_autoastrometry_single(
                img_path=temp_path,
                output_dir=sextractor_out_dir,
                write_crosscheck_files=self.write_crosscheck_files,
                overwrite=True,
                catalog=self.catalog,
                pixel_scale=self.pixel_scale,
                inv=self.inv,
                pa=self.pa,
            )

            # Load up temp path image.header, then delete
            image = self.open_fits(temp_path)

            batch[i] = image

            temp_path.unlink()
            logger.info(
                f"Loaded updated header, and deleted temporary file {temp_path}"
            )

        return batch
