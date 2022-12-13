"""
Central module hosting all shared paths/directory conventions/keys/variables
"""
import logging
import os
import shutil
from collections.abc import Callable
from glob import glob
from pathlib import Path

import toml

logger = logging.getLogger(__name__)

winter_code_dir: Path = Path(__file__).parent.parent.resolve()

toml_path: Path = winter_code_dir.joinpath("pyproject.toml")

with open(toml_path.as_posix(), "r", encoding="UTF-8") as f:
    toml_info = toml.loads(f.read())

package_name: str = toml_info["tool"]["poetry"]["name"]
__version__: str = toml_info["tool"]["poetry"]["version"]

doc_dir = winter_code_dir.joinpath("docs/")

_n_cpu: int | None = os.cpu_count()
if _n_cpu is None:
    default_n_cpu: int = 1
else:
    default_n_cpu = int(_n_cpu / 2)
max_n_cpu: int = int(os.getenv("MAX_N_CPU", max(int(default_n_cpu / 2), 1)))

USE_CACHE: bool = bool(os.getenv("USE_WINTER_CACHE", "true"))

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

raw_img_sub_dir: str = "raw"


def raw_img_dir(
    sub_dir: str = "", raw_dir: Path = base_raw_dir, img_sub_dir: str = raw_img_sub_dir
) -> str:
    """
    Get directory for raw images

    :param sub_dir: sub-dir (night)
    :param raw_dir: Root raw directory for data
    :param img_sub_dir: Default 'raw'
    :return: Full path of raw images
    """
    return os.path.join(raw_dir, os.path.join(str(sub_dir), img_sub_dir))


def get_preprocess_path(raw_img_path: str) -> str:
    return raw_img_path.replace("/raw/", "/preprocess/")


def get_output_dir(
    dir_root: str, sub_dir: str | int = "", output_dir: Path = base_output_dir
) -> Path:
    return output_dir.joinpath(os.path.join(str(sub_dir), dir_root))


def get_output_path(
    base_name: str,
    dir_root: str,
    sub_dir: str | int = "",
    output_dir: Path = base_output_dir,
) -> Path:
    return get_output_dir(
        dir_root, sub_dir=str(sub_dir), output_dir=output_dir
    ).joinpath(base_name)


CACHE_DIR = base_output_dir.joinpath(f"{package_name}_cache")

if not CACHE_DIR.exists():
    CACHE_DIR.mkdir(parents=True)


cal_output_sub_dir = "calibration"


def reduced_img_dir(
    sub_dir: str | int = "", output_dir: Path = base_output_dir
) -> Path:
    return get_output_dir("redux", sub_dir=str(sub_dir), output_dir=output_dir)


def get_mask_path(
    img_path: str | Path,
) -> Path:
    return Path(img_path).with_suffix(".mask.fits")


def get_temp_path(output_dir: Path, file_path: Path) -> Path:
    return output_dir.joinpath("temp_" + file_path.name)


def get_untemp_path(temp_path: Path) -> Path:
    return temp_path.with_name(temp_path.name.replace("temp_", ""))


def copy_temp_file(output_dir: Path, file_path: Path) -> Path:
    output_path = get_temp_path(output_dir=output_dir, file_path=file_path)
    logger.debug(f"Copying from {file_path} to {output_path}")
    shutil.copyfile(file_path, output_path)
    return output_path


raw_img_key = "RAWPATH"
base_name_key = "BASENAME"
ref_img_key = "REFPATH"
proc_history_key = "CALSTEPS"
proc_fail_key = "PROCFAIL"
latest_save_key = "SAVEPATH"
latest_mask_save_key = "MASKPATH"
saturate_key = "SATURATE"
sextractor_header_key = "SRCCAT"
psfex_header_key = "PSFCAT"
norm_psfex_header_key = "NPSFCAT"
ref_psf_key = "REFPSF"
flat_frame_key = "FLATNAME"
bias_frame_key = "BIASNAME"
dark_frame_key = "DARKNAME"
coadd_key = "COADDS"
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
    coadd_key,
    proc_history_key,
    proc_fail_key,
    raw_img_key,
    base_name_key,
]

watchdog_email_key: str = "WATCHDOG_EMAIL"
watchdog_password_key: str = "WATCHDOG_EMAIL_PASSWORD"
watchdog_recipient_key: str = "WATCHDOG_EMAIL_RECIPIENTS"
