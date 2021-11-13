import logging
from winterdrp.telescopes.wirc import update_wirc

logger = logging.getLogger(__name__)

def parse_telescope(header):
    
    telname = None
    
    if header['DETNAM'] == 'Hawaii-II':
        print(header['FPA'])
        if header['FPA'] == 'WIRC':
            telname = "WIRC"
            
    if telname is None:
        err = f"Unable to recognise instrument from fits header. Make sure the instrument has been implemented."
        logger.error(err)
        raise NotImplementedError(err)
        
    logger.debug(f'Telescope recognised as {telname}')
        
    return telname

update_functions = {
    "WIRC": update_wirc,
}

def reformat_raw_data(img):
    
    header = img[0].header
    
    return update_functions[parse_telescope(header)](img)