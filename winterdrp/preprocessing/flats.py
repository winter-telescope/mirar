import numpy as np
import os
from winterdrp.io import open_fits, create_fits
import logging
from glob import glob
from winterdrp.preprocessing.bias import load_master_bias

logger = logging.getLogger(__name__)

base_mflat_name = "master_flat"


def mflat_name(filtername):
    return f"{base_mflat_name}_{filtername}.fits"


def make_master_flats(flatlist, xlolim=500, xuplim=3500, ylolim=500, yuplim=3500, cal_dir='cals'):
    
    if len(flatlist) > 0:

        logger.info(f'Found {len(flatlist)} flat frames')

        with open_fits(flatlist[0]) as img:
            header = img[0].header
            
        nx = header['NAXIS1']
        ny = header['NAXIS2']
        
        master_bias = load_master_bias(cal_dir, header)
        
        filterlist = []
        
        for flat in flatlist:
            with open_fits(flat) as img:
                header = img[0].header
            filterlist.append(header['FILTER'])
            
        flatlist = np.array(flatlist)
            
        for f in list(set(filterlist)):
            
            mask = np.array([x == f for x in flatlist])
            
            cutflatlist = flatlist[mask]
            
            nframes = np.sum(mask)
            
            logger.info(f'Found {nframes} frames for filer {f}')

            flats = np.zeros((ny,nx,nframes))

            for i, flat in enumerate(cutflatlist):
                logger.debug(f'Reading flat {i+1}/{nframes}')

                with open_fits(flat) as img:
                    data = img[0].data
                data = data - master_bias
                median = np.nanmedian(data[xlolim:xuplim,ylolim:yuplim])
                flats[:, :, i] = data/median

            logger.info(f'Median combining {nframes} flats')

            master_flat = np.nanmedian(flats,axis=2)

            with open_fits(flatlist[0]) as img:
                primary_header = img[0].header

            proc_hdu = create_fits(master_flat)  # Create a new HDU with the processed image data
            proc_hdu.header = primary_header       # Copy over the header from the raw file
            proc_hdu.header.add_history('Stacked flat-fielded')

            primary_header['BZERO'] = 0

            mflat_path = os.path.join(cal_dir, mflat_name(f))

            logger.info(f"Saving stacked 'master flat' for filter {f} to {mflat_path}")

            proc_hdu.writeto(mflat_path, overwrite=True)
        
        return 0
      
    else:
        logger.warning("No flat images provided. Proceeding without flat-fielding correction.")


def select_sky_flats():
    raise NotImplementedError


def load_master_flats(cal_dir, header=None):
    
    master_flat_paths = glob(f'{cal_dir}/{base_mflat_name}*.fits')
            
    if len(master_flat_paths) == 0:
        
        try:
            nx = header['NAXIS1']
            ny = header['NAXIS2']

            master_flats = np.zeros((ny,nx))
            
            logger.warning("No master flat found. No flat-fielding correction will be applied.")
        
        except (TypeError, KeyError) as e:
            err = "No master flat files found, and no header info provided to create a dummy image."
            logger.error(err)
            raise FileNotFoundError(err)
            
    else:
        
        master_flats = dict()
        
        for mfpath in master_flat_paths:
                
            f = os.path.basename(mfpath).split("_")[-1][:-5]

            with open_fits(mfpath) as img:
                master_flats[f] = img[0].data
                
    return master_flats


def select_master_flat(all_master_flats, header=None, flat_nan_threshold=0.0):
    
    if isinstance(all_master_flats, np.ndarray):
        master_flat = all_master_flats
    
    elif isinstance(all_master_flats, dict):
        f = header['FILTER']
        try:
            master_flat = all_master_flats[f]
        except KeyError:
            err = f"Unrecognised key {f}. Available filters are {all_master_flats.keys()}"
            logger.error(err)
            raise KeyError(err)
    else:
        err = f"Unrecognised Type for all_master_flats ({type(all_master_flats)}). " \
              f"Was expecting 'numpy.ndarray' or 'dict'."
        logger.error(err)
        raise TypeError(err)
        
     # Mask pixels below a threshold
    
    masked_mflat = np.copy(master_flat)
    if np.any(master_flat < flat_nan_threshold):
        masked_mflat[master_flat < flat_nan_threshold] = np.nan
    
    return masked_mflat