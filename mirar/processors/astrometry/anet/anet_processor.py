"""
Module containing a processor to run astrometry.net
"""
import logging
import os
from pathlib import Path
from typing import Optional

from astropy.io import fits

from mirar.data import ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    get_output_dir,
    get_temp_path,
)
from mirar.processors.astrometry.anet.anet import run_astrometry_net_single
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

ASTROMETRY_TIMEOUT = 900  # astrometry cmd execute timeout, in seconds


class AstrometryNetError(ProcessorError):
    """
    Class for errors in astrometry.net
    """


class AstrometryNet(BaseImageProcessor):
    """Processor to run astrometry.net"""

    base_key = "a-net"

    def __init__(
        self,
        output_sub_dir: str,  # = "a-net"
        scale_bounds: Optional[tuple | list] = None,
        # limits on scale (lower, upper)
        scale_units: Optional[str] = None,  # scale units ('degw', 'amw')
        downsample: Optional[float | int] = None,
        timeout: Optional[
            float
        ] = ASTROMETRY_TIMEOUT,  # astrometry cmd execute timeout, in seconds
        use_sextractor: bool = False,
        sextractor_path: str = "sex",
        search_radius_deg: float = 5.0,
        parity: str = None,
        sextractor_config_path: str = None,
        x_image_key: str = "X_IMAGE",
        y_image_key: str = "Y_IMAGE",
        sort_key_name: str = "FLUX_AUTO",
        use_weight: bool = True,
    ):
        """
        :param output_sub_dir: subdirectory to output astrometry.net results
        :param scale_bounds: limits on scale (lower, upper)
        :param scale_units: scale units ('degw', 'amw')
        :param downsample: downsample by factor of __
        :param timeout: astrometry cmd execute timeout, in seconds
        :param use_sextractor: use sextractor to find sources
        :param sextractor_path: path to sextractor executable (e.g. sex)
        :param search_radius_deg: search radius in degrees
        :param parity: parity of the image, if known (e.g. "odd" or "even")
        :param sextractor_config_path: path to sextractor config file, NOTE that you
        cannot specify other config files (param, conv, nnw, etc.)to astrometry-net.
        Make sure to set the config file to use the correct filter, etc.
        :param x_image_key: key for x-image coordinate in sextractor catalog
        :param y_image_key: key for y-image coordinate in sextractor catalog
        :param sort_key_name: key for sorting sextractor catalog
        """
        super().__init__()

        self.output_sub_dir = output_sub_dir
        self.scale_bounds = scale_bounds
        self.scale_units = scale_units
        self.downsample = downsample
        self.timeout = timeout
        self.use_sextractor = use_sextractor
        self.sextractor_path = sextractor_path
        self.search_radius_deg = search_radius_deg
        self.parity = parity

        self.x_image_key = x_image_key
        self.y_image_key = y_image_key
        self.sort_key_name = sort_key_name
        self.use_weight = use_weight
        self.sextractor_config_path = sextractor_config_path

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
            sextractor_path = f"{self.sextractor_path}"
            if self.use_sextractor & self.use_weight:
                weight_image = image[LATEST_WEIGHT_SAVE_KEY]
                sextractor_path = (
                    f"{self.sextractor_path} -WEIGHT_TYPE MAP_WEIGHT"
                    + f" -WEIGHT_IMAGE {weight_image}"
                )

            run_astrometry_net_single(
                img=temp_path,
                output_dir=anet_out_dir,
                scale_bounds=self.scale_bounds,
                scale_units=self.scale_units,
                downsample=self.downsample,
                timeout=self.timeout,
                use_sextractor=self.use_sextractor,
                sextractor_path=sextractor_path,
                sextractor_config_path=self.sextractor_config_path,
                search_radius_deg=self.search_radius_deg,
                parity=self.parity,
                x_image_key=self.x_image_key,
                y_image_key=self.y_image_key,
                sort_key_name=self.sort_key_name,
            )

            newname = anet_out_dir.joinpath(Path(str(temp_path).split("temp_")[1]))
            if not newname.exists():
                raise AstrometryNetError(
                    f"AstrometryNet did not run successfully - no output "
                    f"file {newname} found."
                )
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
