"""
Module containing a processor to run astrometry.net
"""

import logging
from pathlib import Path
from typing import Callable, Optional

from astropy.io import fits
from astropy.table import Table

from mirar.data import Image, ImageBatch
from mirar.data.utils import write_regions_file
from mirar.errors.exceptions import ProcessorError
from mirar.io import open_fits
from mirar.paths import BASE_NAME_KEY, get_output_dir, get_temp_path
from mirar.processors.astromatic.sextractor.settings import (
    default_conv_path,
    default_starnnw_path,
    parse_sextractor_config,
    write_param_file,
    write_sextractor_config_to_file,
)
from mirar.processors.astrometry.anet.anet import (
    ASTROMETRY_TIMEOUT,
    run_astrometry_net_single,
)
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class AstrometryNetNoSolvedError(ProcessorError):
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
        timeout: float = ASTROMETRY_TIMEOUT,  # astrometry cmd execute timeout, in secs
        use_sextractor: bool = False,
        sextractor_path: str = "sex",
        search_radius_deg: float = 5.0,
        parity: str | None = None,
        sextractor_config_path: str | Path | Callable[[Image], str] | None = None,
        sextractor_params_path: str | Path | None = None,
        sextractor_conv_path: str | Path | None = default_conv_path,
        sextractor_starnnw_path: str | Path | None = default_starnnw_path,
        x_image_key: str = "X_IMAGE",
        y_image_key: str = "Y_IMAGE",
        sort_key_name: str = "MAG_AUTO",
        use_weight: bool = True,
        write_regions: bool = True,
        cache: bool = False,
        no_tweak: bool = False,
    ):
        """
        :param output_sub_dir: subdirectory to output astrometry.net results
        :param scale_bounds: limits on scale (lower, upper)
        :param scale_units: scale units ('degw', 'amw')
        :param downsample: downsample by factor of __
        :param timeout: astrometry cmd execute timeout, in seconds
        :param use_sextractor: use sextractor to find sources
        :param sextractor_path: path to sextractor executable (e.g. sex)
        :param sextractor_params_path: path to sextractor param file
        :param sextractor_conv_path: path to sextractor conv file
        :param sextractor_starnnw_path: path to sextractor starnnw file
        :param search_radius_deg: search radius in degrees
        :param parity: parity of the image, if known (e.g. "odd" or "even")
        :param sextractor_config_path: path to sextractor config file, or a function
        that generates the config file path for a given image. NOTE that you
        cannot specify other config files (param, conv, nnw, etc.)to astrometry-net.
        Make sure to set the config file to use the correct filter, etc.
        :param x_image_key: key for x-image coordinate in sextractor catalog, defaults
        to X_IMAGE, the default from astrometry.net
        :param y_image_key: key for y-image coordinate in sextractor catalog, defaults
        to Y_IMAGE, the default from astrometry.net
        :param sort_key_name: key for sorting sextractor catalog, defaults
        to MAG_AUTO, the default from astrometry.net
        :param no_tweak: do not calculate SIP distortion corrections
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
        if isinstance(sextractor_config_path, str):
            self.sextractor_config_path = Path(sextractor_config_path)
        else:
            self.sextractor_config_path = sextractor_config_path
        self.sextractor_params_path = (
            Path(sextractor_params_path) if sextractor_params_path is not None else None
        )
        self.sextractor_conv_path = (
            Path(sextractor_conv_path) if sextractor_conv_path is not None else None
        )
        self.sextractor_starnnw_path = (
            Path(sextractor_starnnw_path)
            if sextractor_starnnw_path is not None
            else None
        )

        self.write_regions = write_regions

        self.cache = cache
        self.no_tweak = no_tweak

    def description(self) -> str:
        return (
            "Processor to perform astrometric calibration "
            "locally with astrometry.net."
        )

    def get_anet_output_dir(self) -> Path:
        """
        Get the directory to output

        :return: output directory
        """
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    def setup_sextractor_config(self, image: Image) -> (Path | None, list[Path]):
        """
        Setup sextractor config file
        """
        sextractor_temp_files = []
        if isinstance(self.sextractor_config_path, Callable):
            sextractor_config_path = Path(self.sextractor_config_path(image))
        else:
            sextractor_config_path = self.sextractor_config_path
        sextractor_params_path = self.sextractor_params_path
        anet_out_dir = self.get_anet_output_dir()
        if sextractor_config_path is not None:
            if self.sextractor_params_path is None:
                temp_params_path = get_temp_path(
                    anet_out_dir, f"{image[BASE_NAME_KEY]}_sex_anet" f".params"
                )
                sextractor_params_path = temp_params_path
                logger.debug(f"Writing parameters file to " f"{sextractor_params_path}")
                write_param_file(
                    temp_params_path.as_posix(),
                    params=[self.x_image_key, self.y_image_key, self.sort_key_name],
                )
                sextractor_temp_files.append(temp_params_path)
            else:
                assert (
                    sextractor_params_path.exists()
                ), f"sextractor params file {sextractor_params_path} does not exist"
            sextractor_config_dict = parse_sextractor_config(sextractor_config_path)
            sextractor_config_dict["PARAMETERS_NAME"] = sextractor_params_path
            sextractor_config_dict["FILTER_NAME"] = self.sextractor_conv_path
            sextractor_config_dict["STARNNW_NAME"] = self.sextractor_starnnw_path
            temp_config_path = get_temp_path(
                anet_out_dir, f"{image[BASE_NAME_KEY]}_sex_astrom_anet.sex"
            )
            write_sextractor_config_to_file(
                config_dict=sextractor_config_dict,
                config_filename=temp_config_path.as_posix(),
            )
            sextractor_temp_files.append(temp_config_path)
            sextractor_config_path = temp_config_path
        return sextractor_config_path, sextractor_temp_files

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        anet_out_dir = self.get_anet_output_dir()
        anet_out_dir.mkdir(parents=True, exist_ok=True)

        assert len(batch) > 0, "Batch must contain at least one image"

        # Ensure that if a source-extractor config file is provided, it has the
        # correct PARAMETERS_NAME, FILTER_NAME and STARNNW_NAME.
        sextractor_config_path, sextractor_temp_files = self.setup_sextractor_config(
            batch[0]
        )

        for i, image in enumerate(batch):
            base_name = Path(image[BASE_NAME_KEY])
            new_img_path = anet_out_dir.joinpath(base_name)

            temp_path = get_temp_path(anet_out_dir, base_name)
            self.save_fits(image, temp_path, compress=False)

            temp_files = [temp_path, new_img_path]

            sextractor_path = f"{self.sextractor_path}"
            if self.use_sextractor & self.use_weight:

                weight_path = self.save_mask_image(image, temp_path, compress=False)
                temp_files.append(Path(weight_path))

                sextractor_path = (
                    f"{self.sextractor_path} -WEIGHT_TYPE MAP_WEIGHT"
                    + f" -WEIGHT_IMAGE {weight_path}"
                )

            run_astrometry_net_single(
                img_path=temp_path,
                output_dir=anet_out_dir,
                scale_bounds=self.scale_bounds,
                scale_units=self.scale_units,
                downsample=self.downsample,
                timeout=self.timeout,
                use_sextractor=self.use_sextractor,
                sextractor_path=sextractor_path,
                sextractor_config_path=sextractor_config_path,
                search_radius_deg=self.search_radius_deg,
                parity=self.parity,
                x_image_key=self.x_image_key,
                y_image_key=self.y_image_key,
                sort_key_name=self.sort_key_name,
                no_tweak=self.no_tweak,
            )

            anet_output_filees_basepath = new_img_path.as_posix().replace(".fits", "")
            for suffix in [
                ".axy",
                "-objs.png",
                "-ngc.png",
                "-indx.png",
                "-indx.xyls",
                ".corr",
                ".rdls",
                ".match",
                ".solved",
                ".wcs",
            ]:
                temp_files.append(Path(f"{anet_output_filees_basepath}{suffix}"))

            if self.write_regions:
                coords_file = anet_out_dir.joinpath(
                    image[BASE_NAME_KEY].replace(".fits", ".axy")
                )
                if not coords_file.exists():
                    logger.warning(f"Failed to find coords file {coords_file}")
                else:
                    regions_path = anet_out_dir.joinpath(image[BASE_NAME_KEY] + ".reg")
                    logger.debug(f"Loading coords from {coords_file}")
                    with fits.open(coords_file) as hdul:
                        coords_table = Table(hdul[1].data)  # pylint: disable=no-member
                        write_regions_file(
                            regions_path=regions_path,
                            x_coords=coords_table[self.x_image_key],
                            y_coords=coords_table[self.y_image_key],
                            system="image",
                        )

            solved_path = new_img_path.with_suffix(".solved")

            if not solved_path.exists():
                raise AstrometryNetNoSolvedError(
                    f"AstrometryNet did not run successfully - no output "
                    f"file {solved_path} found."
                )

            # Clean up!
            data, hdr = open_fits(new_img_path)
            if "HISTORY" in hdr:
                del hdr["HISTORY"]

            batch[i] = Image(data=data, header=hdr)

            if not self.cache:
                for temp_file in temp_files:
                    temp_file.unlink(missing_ok=True)
                    logger.debug(f"Deleted temporary file {temp_file}")

        if not self.cache:
            for temp_file in sextractor_temp_files:
                temp_file.unlink(missing_ok=True)
                logger.debug(f"Deleted temporary file {temp_file}")

        return batch
