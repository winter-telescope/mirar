import numpy as np
import os
from astropy.io import fits
import logging
from winterdrp.preprocessing.bias import base_mbias_name

logger = logging.getLogger(__name__)

base_mdark_name = "master_dark"

def mdark_name(exptime, norm=False):
    if norm:
        return f"{base_mdark_name}_normed.fits"
    else:
        return f"{base_mdark_name}_{exptime}.fits"

def make_master_dark(darklist, xlolim=500, xuplim=3500, ylolim=500, yuplim=3500, cal_dir='cals', norm=False):
    
    if len(darklist) > 0:

        logger.info(f'Found {len(darklist)} dark frames')

        img = fits.open(darklist[0])
        header = img[0].header
        img.close()
        nx = header['NAXIS1']
        ny = header['NAXIS2']
        exptime = header['EXPTIME']

        darks = np.zeros((ny,nx,len(darklist)))
        
        # Try to load bias image
        
        try:

            img = fits.open(os.path.join(cal_dir, base_mbias_name))
            master_bias = img[0].data
            header = img[0].header
            img.close() 
            
        except FileNotFoundError:
            
            logger.warning("No master bias found. No bias correction will be applied.")
        
            master_bias = np.zeros((ny,nx))
        
        nframes = len(darklist)

        for i, dark in enumerate(darklist):
            logger.debug(f'Reading dark {i+1}/{nframes}')
            img = fits.open(dark)
                        
            if not norm:

                if img[0].header['EXPTIME'] != exptime:
                    logger.error("Error, all dark files do not have the same exposure time."
                                  "Please run with normalise option to create a normalised master dark.")
                    return -1

                darks[:,:,i] = img[0].data - master_bias
                
            else:
                
                dark_exptime = img[0].header['EXPTIME']

                darks[:,:,i] = (img[0].data - master_bias)/dark_exptime
                
            img.close()

        logger.info(f'Median combining {nframes} darks')
            
        master_dark = np.nanmedian(darks,axis=2)

        img = fits.open(darklist[0])
        primaryHeader = img[0].header
        img.close()
        procHDU = fits.PrimaryHDU(master_bias)  # Create a new HDU with the processed image data
        procHDU.header = primaryHeader       # Copy over the header from the raw file
        if norm:
            procHDU.header.add_history('median stacked normalised dark')
            procHDU.header['EXPTIME'] = 1.0

        else:
            procHDU.header.add_history('median stacked dark')
            procHDU.header['EXPTIME'] = exptime

        mdark_path = os.path.join(cal_dir, mdark_name(exptime, norm=norm))
        
        logger.info(f"Saving stacked 'master dark' to {mdark_path}")
        
        procHDU.writeto(mdark_path, overwrite=True)
        return 0
      
    else:
        logger.warning("No dark images provided. Proceeding without dark correction.")
