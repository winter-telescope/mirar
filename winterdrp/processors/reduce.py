import os
import logging
import numpy as np
from winterdrp.io import create_fits
from winterdrp.paths import reduced_img_dir, reduced_img_path, cal_output_dir
from winterdrp.processors.bias import load_master_bias
from winterdrp.processors.dark import load_master_darks, select_master_dark
from winterdrp.processors.flat import load_master_flats, select_master_flat

logger = logging.getLogger(__name__)


def apply_reduction(raw_images, open_fits, sub_dir="", master_bias=None, master_dark=None, master_flat=None,
                    use_norm_dark=False, flat_nan_threshold=0.1, reprocess=True):
    
    cal_dir = cal_output_dir(sub_dir)

    with open_fits(raw_images[0]) as img:
        header = img[0].header

    # Try making output directory, unless it exists
    
    output_dir = reduced_img_dir(sub_dir)
    
    try:
        os.makedirs(output_dir)
    except OSError:
        pass   
    
    # Load up calibration images, if not specified already
    
    if master_bias is None:
    
        master_bias = load_master_bias(cal_dir, header)
        
    if master_dark is None:
        
        master_dark = load_master_darks(cal_dir, header, use_norm=use_norm_dark)
    
    if master_flat is None:
        
        master_flat = load_master_flats(cal_dir, header)
        
    nframes = len(raw_images)
    
    proccessed_list = []
        
    # Loop over science images
        
    for i, raw_img_path in enumerate(raw_images):
        
        img_name = os.path.basename(raw_img_path)

        logger.debug(f"Processing image {i+1}/{nframes} ({img_name})")
        
        output_path = reduced_img_path(img_name, sub_dir=sub_dir)
        
        if np.logical_and(os.path.exists(output_path), reprocess is False):
            logger.debug(f"Skipping image {img_name}, because it has already "
                         f"been processed and 'reprocess' is False.")
            continue
        
        with open_fits(raw_img_path) as img:
            data = img[0].data
            header = img[0].header
        
        if header['OBSTYPE'] not in ['science', "object"]:
            logger.debug(f'Obstype is not science, skipping {raw_img_path}')
            continue
        
        mdark = select_master_dark(master_dark, header)
        mflat = select_master_flat(master_flat, header, flat_nan_threshold=flat_nan_threshold)

        # Master dark????
        
        data_redux = (data - master_bias)/mflat
        proc_hdu = create_fits(data_redux)  # Create a new HDU with the processed image data
        proc_hdu.header = header      # Copy over the header from the raw file
        proc_hdu.header.add_history('Bias corrected and flat-fielded') # Add a note to the header
        proc_hdu.header['BZERO'] = 0
        # Write the reduced frame to disk

        logger.info(f"Saving processed image to {output_path}")
        proccessed_list.append(output_path)
        proc_hdu.writeto(output_path, overwrite=True)

    return proccessed_list
        