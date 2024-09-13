"""
Script containing functions to run astrometry.net locally
"""

import logging
import os
from pathlib import Path
from typing import Optional

from astropy.io import fits

from mirar.errors import ProcessorError
from mirar.utils import ExecutionError, TimeoutExecutionError, execute

logger = logging.getLogger(__name__)


ASTROMETRY_TIMEOUT = 900  # astrometry cmd execute timeout, in seconds


class AstrometryNetExecutionError(ProcessorError):
    """
    Class for errors in astrometry.net
    """


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
    timeout: float = ASTROMETRY_TIMEOUT,  # astrometry cmd execute timeout, in seconds
    use_sextractor: bool = False,
    sextractor_path: str = "sex",
    search_radius_deg: float = 5.0,
    parity: str = None,
    sextractor_config_path: str = None,
    x_image_key: str = "X_IMAGE",
    y_image_key: str = "Y_IMAGE",
    sort_key_name: str = "MAG_AUTO",
    no_tweak: bool = False,
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
    )

    if no_tweak:
        cmd += " --no-tweak "

    if scale_bounds is not None:
        cmd += f" --scale-high {max(scale_bounds)} "
        cmd += f" --scale-low {min(scale_bounds)} "

    if scale_units is not None:
        cmd += f"--scale-units {scale_units} "

    if downsample is not None:
        cmd += f"--downsample {downsample} "

    if use_sextractor:
        cmd += f"--use-source-extractor --source-extractor-path '{sextractor_path}' "

    if sextractor_config_path is not None:
        cmd += f"--source-extractor-config {sextractor_config_path} "

    cmd += f"-X {x_image_key} -Y {y_image_key} -s {sort_key_name} --sort-ascending "

    if parity is not None:
        assert parity in ["pos", "neg"]
        cmd += f"--parity {parity} "

    # cmd with a ra, dec first guess (speeds up solution)
    with fits.open(img_path) as hdul:
        header = hdul[0].header  # pylint: disable=no-member

    ra_req, dec_req = None, None

    if "CRVAL1" in header and "CRVAL2" in header:
        ra_req, dec_req = header["CRVAL1"], header["CRVAL2"]

    elif "RA" in header and "DEC" in header:
        ra_req, dec_req = header["RA"], header["DEC"]  # requested ra, dec

    if ra_req is not None and dec_req is not None:
        try:

            cmd_loc = (
                cmd + f"--ra {ra_req} --dec {dec_req} --radius {search_radius_deg} "
            )  # radius takes on units of ra, dec

            logger.debug(
                f"Running a-net with ra,dec guess and timeout {timeout}. \n"
                f"A-net command:\n {cmd_loc}"
            )

            execute(cmd_loc, output_dir, timeout=timeout)

            assert os.path.isfile(img_path), "Astrometry.net did not solve the image."

        except (ExecutionError, TimeoutExecutionError, KeyError, AssertionError):
            logger.debug("Could not run a-net with ra,dec guess.")

    try:
        logger.debug(f"Running a-net without ra,dec guess.\n" f"A-net command:\n {cmd}")

        execute(cmd, output_dir, timeout=timeout)

        if not os.path.isfile(img_path):
            logger.debug("Second attempt failed.")
            err = "Astrometry.net did not solve the image on second attempt."
            logger.error(err)
            raise AstrometryNetExecutionError(err)

    except (ExecutionError, TimeoutExecutionError) as err:
        raise AstrometryNetExecutionError(err) from err

    return img_path
