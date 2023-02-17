"""
Module containing a processor to run anet.net
"""
import logging
from typing import Optional

from winterdrp.data import ImageBatch
from winterdrp.paths import BASE_NAME_KEY, get_output_dir
from winterdrp.processors.anet.anet import run_astrometry_net_single
from winterdrp.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class AstrometryNet(BaseImageProcessor):
    """Processor to run anet.net"""

    base_key = "a-net"

    def __init__(
        self,
        temp_output_sub_dir: str = "a-net",
        scale_bounds: Optional[tuple | list] = None,  # limits on scale (lower, upper)
        scale_units: Optional[str] = None,  # scale units ('degw', 'amw')
        downsample: Optional[float | int] = None,
    ):
        super().__init__()

        self.temp_output_sub_dir = temp_output_sub_dir
        self.scale_bounds = scale_bounds
        self.scale_units = scale_units
        self.downsample = downsample

    def __str__(self) -> str:
        return "Processor to perform astrometric calibration via anet.net."

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        anet_out_dir = get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

        anet_out_dir.mkdir(parents=True, exist_ok=True)

        for i, image in enumerate(batch):
            temp_path = anet_out_dir.joinpath(image[BASE_NAME_KEY])
            self.save_fits(image, temp_path)

            run_astrometry_net_single(
                img=temp_path,
                output_dir=anet_out_dir,
                scale_bounds=self.scale_bounds,
                scale_units=self.scale_units,
                downsample=self.downsample,
            )

            # Load up temp path image.header, then delete
            image = self.open_fits(temp_path)
            batch[i] = image

            temp_path.unlink()
            logger.info(  # pylint: disable=logging-fstring-interpolation
                f"Loaded updated header, and deleted temporary file {temp_path}"
            )

        return batch
