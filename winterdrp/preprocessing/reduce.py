import os
import logging
import numpy as np
from astropy.io import fits
from winterdrp.paths import raw_img_dir, reduced_img_dir, reduced_img_path, cal_output_dir
from winterdrp.preprocessing.bias import load_master_bias
from winterdrp.preprocessing.dark import load_master_darks, select_master_dark
from winterdrp.preprocessing.flats import load_master_flats, select_master_flat

logger = logging.getLogger(__name__)

def apply_reduction(raw_images, subdir="", master_bias=None, master_dark=None, master_flat=None, use_norm_dark=False, flat_nan_threshold=0.1, reprocess=True):
    
    cal_dir = cal_output_dir(subdir)
    
    with fits.open(raw_images[0]) as img:
        header = img[0].header 
    
    # Try making output directory, unless it exists
    
    output_dir = reduced_img_dir(subdir)
    
    try:
        os.makedirs(output_dir)
    except OSError:
        pass   
    
    # Load up calibration images, if not specified already
    
    if master_bias is None:
    
        master_bias = load_master_bias(cal_dir, header)
        
    if master_dark is None:
        
        all_master_darks = load_master_darks(cal_dir, header, use_norm=use_norm_dark)
    
    if master_flat is None:
        
        all_master_flats = load_master_flats(cal_dir, header)
        
    nframes = len(raw_images)
    
    proccessed_list = []
        
    # Loop over science images
        
    for i, raw_img_path in enumerate(raw_images):
        
        img_name = os.path.basename(raw_img_path)
                
        logger.debug(f"Processing image {i+1}/{nframes} ({img_name})")
        
        output_path = reduced_img_path(subdir, img_name)
        
        if np.logical_and(os.path.exists(output_path), reprocess is False):
            logger.debug(f"Skipping image {img_name}, because it has already been processed and 'reprocess' is False.")
            continue
        
        img = fits.open(raw_img_path)
        data = img[0].data
        header = img[0].header
        
        if header['OBSTYPE'] not in ['science', "object"]:
            logger.debug(f'Obstype is not science, skipping {raw_img_path}')
            continue
        
        master_dark = select_master_dark(all_master_darks, header)
        master_flat = select_master_flat(all_master_flats, header, flat_nan_threshold=flat_nan_threshold)

        # Master dark????
        
        data_redux = (data - master_bias)/master_flat
        #print(procData)
        procHDU = fits.PrimaryHDU(data_redux)  # Create a new HDU with the processed image data
        procHDU.header = header      # Copy over the header from the raw file
        procHDU.header.add_history('Bias corrected and flat-fielded') # Add a note to the header
        primaryHeader['BZERO'] = 0
        # Write the reduced frame to disk

        logging.info(f"Saving processed image to {output_path}")
        proclist.append(output_path)
        procHDU.writeto(output_path, overwrite=True)
        
        