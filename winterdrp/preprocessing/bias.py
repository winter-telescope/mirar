from astropy.io import fits
import os
import numpy as np
import logging

logger = logging.getLogger(__name__)

base_mbias_name = 'master_bias.fits'

def make_master_bias(biaslist, xlolim=200, xuplim=1300, ylolim=200, yuplim=1300, cal_dir='cals'):

    if len(biaslist) > 0:
        
        logger.info(f'Found {len(biaslist)} bias frames')
    
        img = fits.open(biaslist[0])
        header = img[0].header
        img.close()
        nx = header['NAXIS1']
        ny = header['NAXIS2']

        nframes = len(biaslist)

        biases = np.zeros((ny,nx,nframes))
        for i, biasfile in enumerate(biaslist):
            logger.debug(f'Reading bias {i+1}/{nframes}')
            img = fits.open(biasfile)
            biases[:,:,i] = img[0].data
            img.close()

        logger.info(f'Median combining {nframes} biases')

        master_bias = np.nanmedian(biases,axis=2)

        img = fits.open(biaslist[0])
        primaryHeader = img[0].header
        img.close()
        procHDU = fits.PrimaryHDU(master_bias)  # Create a new HDU with the processed image data
        procHDU.header = primaryHeader       # Copy over the header from the raw file
        procHDU.header.add_history('median stacked bias')

        mbias_path = os.path.join(cal_dir, base_mbias_name)

        logger.info(f"Saving stacked 'master bias' to {mbias_path}")

        procHDU.writeto(mbias_path,overwrite=True)
        
    else:
        logger.warning("No bias images provided. Proceeding without bias correction.")
