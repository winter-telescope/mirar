from astropy.io import fits
from winterdrp.telescopes import reformat_raw_data

def read_fits(path):
    with fits.open(path) as raw_img:
        img = reformat_raw_data(raw_img)
    return img

