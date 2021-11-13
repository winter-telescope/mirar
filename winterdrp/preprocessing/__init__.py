import os
import logging
from winterdrp.preprocessing.bias import make_master_bias
from winterdrp.preprocessing.dark import make_master_dark
from winterdrp.preprocessing.flats import make_master_flats, select_sky_flats
from winterdrp.paths import raw_img_dir, cal_output_dir, parse_image_list

logger = logging.getLogger(__name__)

def make_calibration_images(date, flat_method=None):
    
    object_list = parse_image_list(date)

    cal_dir = cal_output_dir(date)
    
    # Make calibration directory, unless it already exists

    try:
        os.makedirs(cal_dir)
    except OSError:
        pass

    make_master_bias(object_list["bias"], cal_dir=cal_dir)
    make_master_dark(object_list["dark"], cal_dir=cal_dir)
    
    if flat_method is None:
        logger.warn("Flat-fielding method not specified. No master flats will be built.")
        
    else:
        if flat_method in ["dome"]:
            flats = object_list["flat"]
        elif flat_method in ["sky"]:
            flats = select_sky_flats(date)
        else:
            msg = f"Selected 'flat_method' ({flat_method}) not recognised. Please specify 'dome' or 'sky, or None to skip this step."
            
            logger.error(msg)
            raise ValueError(msg)

        make_master_flats(flats, cal_dir=cal_dir)