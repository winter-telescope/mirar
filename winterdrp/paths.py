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

def raw_img_dir(date):
    return os.path.join(base_raw_dir, f"{date}/raw")

def cal_output_dir(date):
    return os.path.join(base_output_dir, f"{date}/cals")

def parse_image_list(date):
    
    object_dict = dict()
        
    img_list = glob(f'{raw_img_dir(date)}/*.fits')
    
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
