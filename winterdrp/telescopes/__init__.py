import logging
from winterdrp.telescopes.wirc import update_wirc

logger = logging.getLogger(__name__)


def parse_telescope(header):
    
    telname = None

    try:
        if header['DETNAM'] == 'Hawaii-II':
            if header['FPA'] == 'WIRC':
                telname = "WIRC"
    except KeyError:
        pass
            
    if telname is None:
        msg = f"Unable to recognise instrument from fits header. Make sure the instrument has been implemented. " \
              f"Will now wing it without tidying up the data, and hope for the best."
        logger.warning(msg)
        logger.debug(header)
    else:
        logger.debug(f'Telescope recognised as {telname}')
    return telname


update_functions = {
    "WIRC": update_wirc,
    None: lambda x: x
}


def reformat_raw_data(img):
    
    header = img[0].header
    
    return update_functions[parse_telescope(header)](img)