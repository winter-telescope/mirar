"""
Central module hosting all shared paths/directory conventions/keys/variables
"""
import logging
import os
import shutil
from pathlib import Path

import toml

logger = logging.getLogger(__name__)

winter_code_dir = Path(__file__).parent.parent.resolve()

toml_path = winter_code_dir.joinpath("pyproject.toml")

with open(toml_path.as_posix(), "r", encoding="UTF-8") as f:
    toml_info = toml.loads(f.read())

package_name = toml_info["tool"]["poetry"]["name"]
__version__ = toml_info["tool"]["poetry"]["version"]

doc_dir = winter_code_dir.joinpath("docs/")

_n_cpu = os.cpu_count()
if _n_cpu is None:
    default_n_cpu: int = 1
else:
    default_n_cpu = max(int(_n_cpu / 2), 1)
max_n_cpu: int = int(os.getenv("MAX_N_CPU", default_n_cpu))

# Set up default directories

default_dir = Path.home()

_base_raw_dir: str | None = os.getenv("RAW_DATA_DIR")

if _base_raw_dir is None:
    warning = (
        "No raw data directory specified. "
        "Run 'export RAW_DATA_DIR=/path/to/data' to set. "
        "The raw data directory will need to be specified manually for path function."
        f"The raw directory is being set to {default_dir}."
    )
    logger.warning(warning)
    base_raw_dir: Path = default_dir
else:
    base_raw_dir = Path(_base_raw_dir)

_base_output_dir = os.getenv("OUTPUT_DATA_DIR")

if _base_output_dir is None:
    warning = (
        f"No output data directory specified. "
        f"Run 'export OUTPUT_DATA_DIR=/path/to/data' to set this. "
        f"The output directory is being set to {default_dir}."
    )
    logger.warning(warning)
    base_output_dir = default_dir
else:
    base_output_dir = Path(_base_output_dir)


# Set up special directories
TEMP_DIR = base_output_dir.joinpath(f"{package_name}_temp")
TEMP_DIR.mkdir(exist_ok=True)

RAW_IMG_SUB_DIR = "raw"
CAL_OUTPUT_SUB_DIR = "calibration"


def raw_img_dir(
    sub_dir: str = "", raw_dir: Path = base_raw_dir, img_sub_dir: str = RAW_IMG_SUB_DIR
) -> Path:
    """
    Get directory for raw images

    :param sub_dir: sub-dir (night)
    :param raw_dir: Root raw directory for data
    :param img_sub_dir: Default 'raw'
    :return: Full path of raw images
    """
    return raw_dir.joinpath(os.path.join(str(sub_dir), img_sub_dir))


def get_output_dir(
    dir_root: str, sub_dir: str | int = "", output_dir: Path = base_output_dir
) -> Path:
    """
    Generic function to get a full output directory combining dir_root, sub_dir and
    the parent output directory

    :param dir_root: directory within subdir, e.g 'raw' or 'processed'
    :param sub_dir: subdirectory in parent directory, typically a night e.g 20221223
    :param output_dir: parent output directory
    :return: full output directory
    """
    return output_dir.joinpath(os.path.join(str(sub_dir), dir_root))


def get_output_path(
    base_name: str,
    dir_root: str,
    sub_dir: str | int = "",
    output_dir: Path = base_output_dir,
) -> Path:
    """
    Generic function to get a full output path combining the file name, dir_root,
    sub_dir and the parent output directory

    :param base_name: name of file
    :param dir_root: directory within subdir, e.g 'raw' or 'processed'
    :param sub_dir: subdirectory in parent directory, typically a night e.g 20221223
    :param output_dir: parent output directory
    :return: full output directory
    """
    return get_output_dir(
        dir_root, sub_dir=str(sub_dir), output_dir=output_dir
    ).joinpath(base_name)


def get_weight_path(
    img_path: str | Path,
) -> Path:
    """
    Returns a weight image path

    :param img_path: parent image
    :return: custom path for weight image
    """
    return Path(img_path).with_suffix(".weight.fits")


def get_temp_path(output_dir: Path, file_path: Path | str) -> Path:
    """
    Gets a temporary path, in output dir, with name of file_path

    :param output_dir: Output directory
    :param file_path: current path of file
    :return: temporary path
    """
    return output_dir.joinpath("temp_" + Path(file_path).name)


def get_untemp_path(temp_path: Path | str) -> Path:
    """
    Converts a temporary path to a regular path.

    Essentially undoes ..:func:`winterdrp.path.get_temp_path`

    :param temp_path: temporary file path
    :return: normal file path
    """
    temp_path = Path(temp_path)
    return temp_path.with_name(temp_path.name.replace("temp_", ""))


def copy_temp_file(output_dir: Path, file_path: Path) -> Path:
    """
    Copies a file at file_path to a temporary path in output dir,
    then returns temp path

    :param output_dir: output directory
    :param file_path: file to cope
    :return: path of temporary file
    """
    output_path = get_temp_path(output_dir=output_dir, file_path=file_path)
    logger.debug(f"Copying from {file_path} to {output_path}")
    shutil.copyfile(file_path, output_path)
    return output_path

RAW_IMG_KEY = "RAWPATH"
BASE_NAME_KEY = "BASENAME"
REF_IMG_KEY = "REFPATH"
PROC_HISTORY_KEY = "CALSTEPS"
PROC_FAIL_KEY = "PROCFAIL"
LATEST_SAVE_KEY = "SAVEPATH"
LATEST_WEIGHT_SAVE_KEY = "WGHTPATH"
SEXTRACTOR_HEADER_KEY = "SRCCAT"
PSFEX_CAT_KEY = "PSFCAT"
NORM_PSFEX_KEY = "NPSFCAT"
REF_PSF_KEY = "REFPSF"
FLAT_FRAME_KEY = "FLATNAME"
BIAS_FRAME_KEY = "BIASNAME"
DARK_FRAME_KEY = "DARKNAME"
COADD_KEY = "COADDS"
GAIN_KEY = "GAIN"
sextractor_checkimg_keys = {
    "BACKGROUND": "BKGPT",
    "BACKGROUND_RMS": "BKGRMS",
    "MINIBACKGROUND": "MINIBKG",
    "MINIBACK_RMS": "MINIBGRM",
}

core_fields = [
    "OBSCLASS",
    "TARGET",
    "UTCTIME",
    COADD_KEY,
    PROC_HISTORY_KEY,
    PROC_FAIL_KEY,
    RAW_IMG_KEY,
    BASE_NAME_KEY,
]

MONITOR_EMAIL_KEY = "WATCHDOG_EMAIL"
MONITOR_PASSWORD_KEY = "WATCHDOG_EMAIL_PASSWORD"
MONITOR_RECIPIENT_KEY = "WATCHDOG_EMAIL_RECIPIENTS"

