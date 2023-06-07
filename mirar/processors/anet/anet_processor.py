"""
Module containing a processor to run astrometry.net
"""
import logging
import os
from pathlib import Path
from typing import Optional

from astropy.io import fits

from mirar.data import ImageBatch
from mirar.paths import BASE_NAME_KEY, get_output_dir, get_temp_path
from mirar.processors.anet.anet import run_astrometry_net_single
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


ASTROMETRY_TIMEOUT = 900  # astrometry cmd execute timeout, in seconds


class AstrometryNet(BaseImageProcessor):
    """Processor to run astrometry.net"""

    base_key = "a-net"

    def __init__(
        self,
        output_sub_dir: str,  # = "a-net"
        scale_bounds: Optional[tuple | list] = None,  # limits on scale (lower, upper)
        scale_units: Optional[str] = None,  # scale units ('degw', 'amw')
        downsample: Optional[float | int] = None,
        timeout: Optional[
            float
        ] = ASTROMETRY_TIMEOUT,  # astrometry cmd execute timeout, in seconds
    ):
        super().__init__()

        self.output_sub_dir = output_sub_dir
        self.scale_bounds = scale_bounds
        self.scale_units = scale_units
        self.downsample = downsample
        self.timeout = timeout

    def __str__(self) -> str:
        return "Processor to perform astrometric calibration via astrometry.net."

    def get_anet_output_dir(self) -> Path:
        """
        Get the directory to output

        :return: output directory
        """
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        anet_out_dir = self.get_anet_output_dir()
        cache = False
        try:
            os.makedirs(anet_out_dir)
        except OSError:
            pass

        for i, image in enumerate(batch):
            temp_path = get_temp_path(anet_out_dir, image[BASE_NAME_KEY])
            if not os.path.exists(temp_path):
                self.save_fits(image, temp_path)

            temp_files = [temp_path]

            run_astrometry_net_single(
                img=temp_path,
                output_dir=anet_out_dir,
                scale_bounds=self.scale_bounds,
                scale_units=self.scale_units,
                downsample=self.downsample,
                timeout=self.timeout,
            )

            newname = anet_out_dir.joinpath(Path(str(temp_path).split("temp_")[1]))
            solved = fits.open(newname)
            hdr = solved[0].header  # pylint: disable=no-member

            del hdr["HISTORY"]

            fits.writeto(  # pylint: disable=no-member
                newname,
                fits.open(temp_files[0])[0].data,  # pylint: disable=no-member
                hdr,
                overwrite=True,
            )  # pylint: disable=no-member
            batch[i] = self.open_fits(newname)  # pylint: disable=no-member

            if not cache:
                for temp_file in temp_files:
                    os.remove(temp_file)
                    logger.info(f"Deleted temporary file {temp_file}")

        return batch
