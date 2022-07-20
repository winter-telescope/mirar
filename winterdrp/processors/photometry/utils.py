from astropy.io import fits
import numpy as np
import logging

logger = logging.getLogger(__name__)


def make_cutouts(
        image_path: str,
        position,
        half_size
):
    data = fits.getdata(image_path)
    y_image_size, x_image_size = np.shape(data)
    x, y = position
    # logger.debug(f'{x},{y},{np.shape(data)}')
    if np.logical_and(x < half_size,y < half_size):
        cutout = data[0:y+half_size+1, 0:x+half_size+1]
        n_xpix = half_size-y
        n_ypix = half_size-x
        cutout = np.pad(cutout, ((n_ypix, 0), (n_xpix, 0)), 'constant')

    elif np.logical_and(x+half_size+1 > x_image_size, y+half_size+1 > y_image_size):
        cutout = data[y - half_size: y_image_size, x-half_size, x_image_size]
        n_xpix = (half_size+x+1) - x_image_size
        n_ypix = (half_size+y+1) - y_image_size
        cutout = np.pad(cutout, ((0, n_ypix), (0, n_xpix)), 'constant')

    elif y < half_size:
        logger.info(f'Cutout parameters are {y + half_size + 1}, {x - half_size}, {x + half_size + 1},{y_image_size},'
                    f'{x_image_size}')
        cutout = data[0:y + half_size + 1, x - half_size:x + half_size + 1]
        n_pix = half_size - y
        cutout = np.pad(cutout, ((n_pix, 0), (0, 0)), 'constant')

    elif y + half_size + 1 > y_image_size:
        cutout = data[y - half_size: y_image_size, x - half_size: x + half_size + 1]
        n_pix = (half_size + y + 1) - y_image_size
        cutout = np.pad(cutout, ((0, n_pix), (0, 0)), 'constant')

    elif x < half_size:
        cutout = data[y - half_size: y + half_size + 1, 0:x + half_size + 1]
        n_pix = half_size - x
        cutout = np.pad(cutout, ((0, 0), (n_pix, 0)), 'constant')
    elif x + half_size > x_image_size:
        cutout = data[y - half_size:y + half_size + 1, x - half_size:x_image_size]
        n_pix = (half_size + x + 1) - x_image_size
        cutout = np.pad(cutout, ((0, 0), (0, n_pix)), 'constant')
    else:
        cutout = data[y - half_size:y + half_size + 1, x - half_size:x + half_size + 1]
    return cutout
