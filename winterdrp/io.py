from astropy.io import fits
from winterdrp.telescopes import reformat_raw_data


def open_fits(path):
    img = fits.open(path)
    img = reformat_raw_data(img)
    return img


def create_fits(data):
    return fits.PrimaryHDU(data)
