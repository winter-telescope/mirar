import os
from astropy.io import fits
from glob import glob
import logging
from collections.abc import Callable

logger = logging.getLogger(__name__)

base_raw_dir = os.getenv("RAW_DATA_DIR")

if base_raw_dir is None:
    err = "No raw data directory specified. Run 'export RAW_DATA_DIR=/path/to/data'"
    logger.error(err)
    raise ValueError(err)
    
base_output_dir = os.getenv("OUTPUT_DATA_DIR")

if base_output_dir is None:
    err = "No output data directory specified. Run 'export OUTPUT_DATA_DIR=/path/to/data'"
    logger.error(err)
    raise ValueError(err)


def raw_img_dir(
        sub_dir=""
):
    return os.path.join(base_raw_dir, os.path.join(sub_dir, "raw"))


def output_dir(
        dir_root: str,
        sub_dir: str | int = ""
):
    return os.path.join(base_output_dir, os.path.join(sub_dir, dir_root))


def output_path(
        base_name: str,
        dir_root: str,
        sub_dir: str | int = ""
):
    return os.path.join(output_dir(dir_root, sub_dir=sub_dir), base_name)


def cal_output_dir(
        sub_dir: str | int = ""
):
    return output_dir("cals", sub_dir=sub_dir)


def reduced_img_dir(
        sub_dir: str | int = ""
):
    return output_dir("redux", sub_dir=sub_dir)


def reduced_img_path(
        img_name: str,
        sub_dir: str | int = ""
):
    return os.path.join(reduced_img_dir(sub_dir), img_name)


def observing_log_dir(
        sub_dir: str | int = ""
):
    return os.path.join(base_raw_dir, sub_dir)


def astrometry_output_dir(
        sub_dir: str = "",
        astro_pass: int = 1
):
    return output_dir(f"astrometry_{astro_pass}", sub_dir=sub_dir)


def parse_image_list(
        sub_dir: str | int = "",
        group_by_object: bool = True,
        base_dir_f: Callable[[str], str] = raw_img_dir
):
    
    object_dict = dict()
        
    img_list = glob(f'{base_dir_f(sub_dir)}/*.fits')

    if not group_by_object:
        return sorted(img_list)
    
    for img_file in img_list:
        img = fits.open(img_file)
        
        obj = img[0].header['OBJECT']
        
        if obj not in object_dict.keys():
            object_dict[obj] = [img_file]
        else:
            object_dict[obj].append(img_file)
            
    logger.debug(f'Data contains {len(object_dict.keys())} objects: {list(object_dict.keys())}')
    
    for key in ["dark", "flat", "bias"]:
        if key not in object_dict.keys():
            object_dict[key] = []

    return object_dict
