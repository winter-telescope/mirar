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

ANET_INDEX_DIR_KEY = "ANET_INDEX_DIR"

DEFAULT_INDEX_SERIES = "4200"

ANET_DOWNLOAD_URL = "http://data.astrometry.net"

# Each series is served from its own directory, apart from 5200 and its lite version
INDEX_SERIES_URLS = {
    "4100": f"{ANET_DOWNLOAD_URL}/4100/",
    "4200": f"{ANET_DOWNLOAD_URL}/4200/",
    "5000": f"{ANET_DOWNLOAD_URL}/5000/",
    "5200": "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/",
    "6000": f"{ANET_DOWNLOAD_URL}/6000/",
    "6100": f"{ANET_DOWNLOAD_URL}/6100/",
}


class AstrometryNetExecutionError(ProcessorError):
    """
    Class for errors in astrometry.net
    """


class AstrometryNetIndexError(ProcessorError):
    """
    Class for missing astrometry.net index files
    """


def get_index_file_pattern(index_series: str) -> str:
    """
    Get the filename pattern of the index files belonging to a series

    A series is named after its catalogue, e.g. '5200', while each file within it is
    named after the angular scale it covers, e.g. 'index-5204-11.fits'

    :param index_series: Series of index files, e.g. '5200'
    :return: Pattern matching the index files of that series
    """
    return f"index-{index_series[:2]}*.fits"


def get_index_dir(index_series: str) -> Path:
    """
    Get the directory holding one series of astrometry.net index files

    The ANET_INDEX_DIR environment variable sets the base directory, and each
    series lives in its own subdirectory of it, mirroring the layout of the
    astrometry.net download server. That keeps the series separate, which
    matters because the solution inherits the reference frame of whichever
    index solved it.

    :param index_series: Series of index files to use, e.g. '5200'
    :return: Directory holding that series of index files
    """
    base_dir = os.getenv(ANET_INDEX_DIR_KEY)

    if base_dir is None:
        err = (
            f"astrometry.net needs index files, which are not included in its "
            f"installation. Set the '{ANET_INDEX_DIR_KEY}' environment variable to "
            f"the directory holding them, with one subdirectory per series, so that "
            f"the '{index_series}' series sits in "
            f"'${ANET_INDEX_DIR_KEY}/{index_series}'. Download it from "
            f"{INDEX_SERIES_URLS.get(index_series, ANET_DOWNLOAD_URL)}"
        )
        logger.error(err)
        raise AstrometryNetIndexError(err)

    index_dir = Path(base_dir) / index_series

    if not index_dir.is_dir():
        raise AstrometryNetIndexError(
            missing_index_error(index_series, index_dir, "There is no such directory")
        )

    pattern = get_index_file_pattern(index_series)

    if len(list(index_dir.glob(pattern))) == 0:
        raise AstrometryNetIndexError(
            missing_index_error(
                index_series, index_dir, f"That directory holds no '{pattern}' files"
            )
        )

    return index_dir


def missing_index_error(index_series: str, index_dir: Path, reason: str) -> str:
    """
    Build the error message for a series of index files that cannot be found

    :param index_series: Series of index files that was requested, e.g. '5200'
    :param index_dir: Directory the index files were expected in
    :param reason: Why the series could not be used, e.g. 'There is no such directory'
    :return: Error message
    """
    err = (
        f"Found no astrometry.net index files for the '{index_series}' series. "
        f"Expected files matching '{get_index_file_pattern(index_series)}' in "
        f"'{index_dir}'. {reason}. The '{ANET_INDEX_DIR_KEY}' environment variable "
        f"sets the base directory, and each series belongs in its own subdirectory "
        f"of it, so download the '{index_series}' series from "
        f"{INDEX_SERIES_URLS.get(index_series, ANET_DOWNLOAD_URL)} and unpack it "
        f"into '{index_dir}'"
    )
    logger.error(err)
    return err


def run_astrometry_net(images: Path | list, output_dir: Path, *args, **kwargs):
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
    img_path: Path,
    output_dir: Path,
    scale_bounds: Optional[tuple | list] = None,  # limits on scale (lower, upper)
    scale_units: Optional[str] = None,  # scale units ('degw', 'amw')
    downsample: Optional[float | int] = None,  # downsample by factor of __
    timeout: float = ASTROMETRY_TIMEOUT,  # astrometry cmd execute timeout, in seconds
    use_sextractor: bool = False,
    sextractor_path: str = "sex",
    search_radius_deg: float = 5.0,
    parity: str | None = None,
    sextractor_config_path: str | None = None,
    x_image_key: str = "X_IMAGE",
    y_image_key: str = "Y_IMAGE",
    sort_key_name: str = "MAG_AUTO",
    no_tweak: bool = False,
    index_series: str = DEFAULT_INDEX_SERIES,
):
    """
    function to run astrometry.net locally on one image, with options to adjust settings
    default: solve-field <img> -D <output_dir> -N <newname> -O
    """
    # name for new file if a-net solves (otherwise a-net writes to '<img>.new')
    newname = output_dir / Path(str(img_path).split("temp_")[1])
    basename = (str(img_path).split("temp_")[1]).split(".fits")[0]

    # run a-net (solve-field)
    cmd = (
        f"solve-field {img_path} "
        f"--dir {output_dir} "
        f"--new-fits {newname} "
        f"--overwrite "
        f"--out {basename} "  # use this base name for outputs (instead of 'temp_...')
        f"--no-plots "
    )

    cmd += f"--index-dir {get_index_dir(index_series)} "

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

    solved = False

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

            solved = True

        except (ExecutionError, TimeoutExecutionError, KeyError, AssertionError):
            logger.debug("Could not run a-net with ra,dec guess.")

    if not solved:
        try:
            logger.debug(
                f"Running a-net without ra,dec guess.\n" f"A-net command:\n {cmd}"
            )

            execute(cmd, output_dir, timeout=timeout)

            if not os.path.isfile(img_path):
                logger.debug("Second attempt failed.")
                err = "Astrometry.net did not solve the image on second attempt."
                logger.error(err)
                raise AstrometryNetExecutionError(err)

        except (ExecutionError, TimeoutExecutionError) as err:
            raise AstrometryNetExecutionError(err) from err

    return img_path
