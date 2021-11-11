import os
from winterdrp.preprocessing.bias import make_master_bias
from winterdrp.preprocessing.dark import make_master_dark
from winterdrp.preprocessing.flats import make_master_flats
from winterdrp.paths import raw_img_dir, cal_output_dir, parse_image_list

def make_calibration_images(date, norm_dark=True):
    
    object_list = parse_image_list(date)

    cal_dir = cal_output_dir(date)
    
    # Make calibration directory, unless it already exists

    try:
        os.makedirs(cal_dir)
    except OSError:
        pass

    make_master_bias(object_list["bias"], cal_dir=cal_dir)
    make_master_dark(object_list["dark"], cal_dir=cal_dir, norm=norm_dark)
    make_master_flats(object_list["flat"], cal_dir=cal_dir)