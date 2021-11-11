import numpy as np
import os
from astropy.io import fits
import logging
from winterdrp.preprocessing.bias import base_mbias_name

logger = logging.getLogger(__name__)

base_mflat_name = "master_flat"

def mflat_name(filtername):
    return f"{base_mflat_name}_{filtername}.fits"

def make_master_flats(flatlist, xlolim=500, xuplim=3500, ylolim=500, yuplim=3500, cal_dir='cals'):
    
    if len(flatlist) > 0:

        logger.info(f'Found {len(flatlist)} flat frames')

        img = fits.open(flatlist[0])
        header = img[0].header
        img.close()
        nx = header['NAXIS1']
        ny = header['NAXIS2']
        
        # Try to load bias image
        
        try:

            img = fits.open(os.path.join(cal_dir, base_mbias_name))
            master_bias = img[0].data
            header = img[0].header
            img.close() 
            
        except FileNotFoundError:
            
            logger.warning("No master bias found. No bias correction will be applied.")
        
            master_bias = np.zeros((ny,nx))
        
        filterlist = []
        
        for flat in flatlist:
            img = fits.open(flat)
            data = img[0].data
            header = img[0].header
            
            filterlist.append(header['FILTER'])
            
        flatlist = np.array(flatlist)
            
        for f in list(set(filterlist)):
            
            mask = flatlist == f
            
            cutflatlist = flatlist[mask]
            
            nframes = np.sum(mask)
            
            logger.info(f'Found {nframes} frames for filer {f}')

            flats = np.zeros((ny,nx,nframes))

            for i, flat in enumerate(cutflatlist):
                logger.debug(f'Reading flat {i+1}/{nframes}')
                img = fits.open(flat)
                data = img[0].data
                img.close()
                data = data - master_bias
                median = np.nanmedian(data[xlolim:xuplim,ylolim:yuplim])
                flats[:,:,i] = data/median

            logger.info(f'Median combining {nframes} flats')

            master_flat = np.nanmedian(flats,axis=2)

            img = fits.open(flatlist[0])
            primaryHeader = img[0].header
            img.close()
            procHDU = fits.PrimaryHDU(masterFlat)  # Create a new HDU with the processed image data
            procHDU.header = primaryHeader       # Copy over the header from the raw file
            procHDU.header.add_history('Stacked flat-fielded')

            primaryHeader['BZERO'] = 0

            mflat_path = os.path.join(cal_dir, mflat_name(f))

            logger.info(f"Saving stacked 'master flat' for filter {f} to {mflat_path}")

            procHDU.writeto(mflat_path, overwrite=True)
        
        return 0
      
    else:
        logger.warning("No flat images provided. Proceeding without flat-fielding correction.")