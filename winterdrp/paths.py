import os
import logging
from astropy.io import fits
from glob import glob
import logging

logger = logging.getLogger(__name__)

base_raw_dir = os.getenv("RAW_DATA_DIR")

if base_raw_dir is None:
    logger.error("No raw data directory specified. Run 'export RAW_DATA_DIR=/path/to/data'")
    raise ValueError("No raw data directory specified.")
    
base_output_dir = os.getenv("OUTPUT_DATA_DIR")

if base_output_dir is None:
    logger.error("No output data directory specified. Run 'export OUTPUT_DATA_DIR=/path/to/data'")
    raise ValueError("No output data directory specified.")

def raw_img_dir(subdir=""):
    return os.path.join(base_raw_dir, f"{subdir}/raw")

def cal_output_dir(subdir=""):
    return os.path.join(base_output_dir, f"{subdir}/cals")

def reduced_img_dir(subdir=""):
    return os.path.join(base_output_dir, f"{subdir}/redux")

def reduced_img_path(img_name, subdir=""):
    return os.path.join(reduced_img_dir(subdir), img_name)
    
def parse_image_list(subdir=""):
    
    object_dict = dict()
        
    img_list = glob(f'{raw_img_dir(subdir)}/*.fits')
    
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
