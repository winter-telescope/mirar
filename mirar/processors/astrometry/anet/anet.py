"""
Script containing functions to run astrometry.net locally
"""
import logging
import os
from pathlib import Path
from typing import Optional

from astropy.io import fits

from mirar.utils import ExecutionError, execute

logger = logging.getLogger(__name__)


class AstrometryNetError(ExecutionError):
    """Astrometry.net-related error"""


def run_astrometry_net(images: str | list, output_dir: str, *args, **kwargs):
    """
    function to execute `run_astrometry_net_single` on several images in batch
    """
    if not isinstance(images, list):
        images = [images]

    # make output directory if it doesn't exist

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    for img in images:
        run_astrometry_net_single(img, output_dir, *args, **kwargs)


def run_astrometry_net_single(
    img_path: str | Path,
    output_dir: str | Path,
    scale_bounds: Optional[tuple | list] = None,  # limits on scale (lower, upper)
    scale_units: Optional[str] = None,  # scale units ('degw', 'amw')
    downsample: Optional[float | int] = None,  # downsample by factor of __
    timeout: Optional[float] = None,  # astrometry cmd execute timeout, in seconds
    use_sextractor: bool = False,
    sextractor_path: str = "sex",
    search_radius_deg: float = 5,
    parity: str = None,
    sextractor_config_path: str = None,
    x_image_key: str = "X_IMAGE",
    y_image_key: str = "Y_IMAGE",
    sort_key_name: str = "MAG_AUTO",
):
    """
    function to run astrometry.net locally on one image, with options to adjust settings
    default: solve-field <img> -D <output_dir> -N <newname> -O
    """
    # name for new file if a-net solves (otherwise a-net writes to '<img>.new')
    newname = output_dir.joinpath(Path(str(img_path).split("temp_")[1]))
    basename = (str(img_path).split("temp_")[1]).split(".fits")[0]

    # run a-net (solve-field)
    cmd = (
        f"solve-field {img_path} "
        f"--dir {output_dir} "
        f"--new-fits {newname} "
        f"--overwrite "
        f"--out {basename} "  # use this base name for outputs (instead of 'temp_...')
        f"--use-source-extractor "
    )

    if scale_bounds is not None:
        cmd += f" --scale-high {max(scale_bounds)} "
        cmd += f" --scale-low {min(scale_bounds)} "

    if scale_units is not None:
        cmd += f"--scale-units {scale_units} "

    if downsample is not None:
        cmd += f"--downsample {downsample} "

    # cmd with a ra, dec first guess (speeds up solution)
    header = fits.open(img_path)[0].header  # pylint: disable=no-member
    ra_req, dec_req = header["RA"], header["DEC"]  # requested ra, dec
    if use_sextractor:
        cmd += f"--use-source-extractor --source-extractor-path '{sextractor_path}' "

    if sextractor_config_path is not None:
        cmd += f"--source-extractor-config {sextractor_config_path} "

    cmd += f"-X {x_image_key} -Y {y_image_key} -s {sort_key_name} --sort-ascending "

    if parity is not None:
        assert parity in ["pos", "neg"]
        cmd += f"--parity {parity} "

    cmd_loc = (
        cmd + f"--ra {ra_req} --dec {dec_req} --radius {search_radius_deg} "
    )  # radius takes on units of ra, dec

    try:
        logger.debug(
            f"Running a-net with ra,dec guess and timeout {timeout}. \n"
            f"A-net command:\n {cmd_loc}"
        )

        if timeout is not None:
            execute(cmd_loc, output_dir, timeout=timeout)
        else:
            execute(cmd_loc, output_dir)

        if not os.path.isfile(img_path):
            logger.debug(
                f"First attempt failed. "
                f"Rerunning a-net without ra,dec guess.\n"
                f"A-net command:\n {cmd}"
            )

            if timeout is not None:
                execute(cmd, output_dir, timeout=timeout)
            else:
                execute(cmd, output_dir)

            if not os.path.isfile(img_path):
                logger.debug("Second attempt failed.")
    except ExecutionError as err:
        raise AstrometryNetError(err) from err

    return img_path
