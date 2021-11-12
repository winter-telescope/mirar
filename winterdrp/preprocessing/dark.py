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
        return f"{base_mdark_name}_{exptime:.0f}s.fits"

def make_master_dark(darklist, xlolim=500, xuplim=3500, ylolim=500, yuplim=3500, cal_dir='cals', norm=False):
    
    if len(darklist) > 0:

        logger.info(f'Found {len(darklist)} dark frames')

        img = fits.open(darklist[0])
        header = img[0].header
        img.close()
        nx = header['NAXIS1']
        ny = header['NAXIS2']

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
        
        # If normed, loop over all darks

        if norm:
            logger.info("Norm is True. Making one 'master_dark' combining all exposure times.")
            dark_loop = [(darklist, 1.0, "Median stacked normalised dark")]
            
        # If not normed, split darks into blocks by exposure time
        
        else:
            logger.info("Norm is False. Making one 'master_dark' for each exposure time.")
            
            explist = []
        
            for dark in darklist:
                img = fits.open(dark)
                data = img[0].data
                header = img[0].header

                explist.append(header['EXPTIME'])
                
            exps = sorted(list(set(explist)))
            
            logger.info(f'Found {len(exps)} different exposure times: {exps}')

            darklist = np.array(darklist)
            
            dark_loop = []

            for exp in exps:

                mask = np.array([x == exp for x in explist])
                
                dark_loop.append((darklist[mask], exp, "Median stacked dark"))
                
        # Loop over each set of darks and create a master_dark
                
        for (cutdarklist, exptime, history) in dark_loop:
        
            nframes = len(cutdarklist)

            for i, dark in enumerate(cutdarklist):

                img = fits.open(dark)
                dark_exptime = img[0].header['EXPTIME']

                logger.debug(f'Reading dark {i+1}/{nframes} with exposure time {dark_exptime}')

                darks[:,:,i] = (img[0].data - master_bias) * exptime/dark_exptime

                img.close()

            master_dark = np.nanmedian(darks,axis=2)

            img = fits.open(darklist[0])
            primaryHeader = img[0].header
            img.close()
            procHDU = fits.PrimaryHDU(master_bias)  # Create a new HDU with the processed image data
            procHDU.header = primaryHeader       # Copy over the header from the raw file

            procHDU.header.add_history(history)
            procHDU.header['EXPTIME'] = exptime

            mdark_path = os.path.join(cal_dir, mdark_name(exptime, norm=norm))

            logger.info(f"Saving stacked 'master dark' combining {nframes} exposures to {mdark_path}")

            procHDU.writeto(mdark_path, overwrite=True)
        return 0
      
    else:
        logger.warning("No dark images provided. Proceeding without dark correction.")
