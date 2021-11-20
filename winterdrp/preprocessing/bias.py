from winterdrp.io import open_fits, create_fits
import os
import numpy as np
import logging

logger = logging.getLogger(__name__)

base_mbias_name = 'master_bias.fits'


def make_master_bias(biaslist, cal_dir='cals'):

    if len(biaslist) > 0:
        
        logger.info(f'Found {len(biaslist)} bias frames')
    
        with open_fits(biaslist[0]) as img:
            header = img[0].header
        
        nx = header['NAXIS1']
        ny = header['NAXIS2']

        nframes = len(biaslist)

        biases = np.zeros((ny, nx, nframes))
        
        for i, bias in enumerate(biaslist):
            logger.debug(f'Reading bias {i+1}/{nframes}')
            with open_fits(bias) as img:
                biases[:, :, i] = img[0].data

        logger.info(f'Median combining {nframes} biases')

        master_bias = np.nanmedian(biases, axis=2)

        with open_fits(biaslist[0]) as img:
            primary_header = img[0].header

        proc_hdu = create_fits(master_bias)  # Create a new HDU with the processed image data
        proc_hdu.header = primary_header      # Copy over the header from the raw file
        proc_hdu.header.add_history('median stacked bias')

        mbias_path = os.path.join(cal_dir, base_mbias_name)

        logger.info(f"Saving stacked 'master bias' to {mbias_path}")

        proc_hdu.writeto(mbias_path, overwrite=True)
        
    else:
        logger.warning("No bias images provided. No master bias created.")


def load_master_bias(cal_dir, header=None):
    
    # Try to load bias image
        
    try:
        with open_fits(os.path.join(cal_dir, base_mbias_name)) as img:
            master_bias = img[0].data

    except FileNotFoundError:
        
        try:
            
            nx = header['NAXIS1']
            ny = header['NAXIS2']

            master_bias = np.zeros((ny, nx))
            
            logger.warning("No master bias found. No bias correction will be applied.")
            
        except KeyError:
            
            err = "No master bias files found, and no header info provided to create a dummy image."
            logger.error(err)
            raise FileNotFoundError(err)

    return master_bias
