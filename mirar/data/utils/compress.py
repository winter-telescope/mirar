"""
Module for encoding/decoding numpy 2D images to gziped fits files
"""

import gzip
import io

import numpy as np
from astropy.io import fits


def decode_img(compressed_bytes: bytes) -> np.ndarray:
    """
    Function to parse a cutout (gziped fits file) into a numpy array

    :param compressed_bytes: Gziped fits file bytes
    :return: Numpy array of the image
    """
    with gzip.open(io.BytesIO(compressed_bytes), "r") as gzipped_f:
        with fits.open(io.BytesIO(gzipped_f.read()), ignore_missing_simple=True) as hdu:
            data = hdu[0].data  # pylint: disable=no-member
    return data


def encode_img(image: np.array):
    """
    make bit images of the cutouts for the marshal
    Args:
        image: 2D numpy array

    Returns:
        compressed bytes
    """
    with io.BytesIO() as buffer_1:
        with io.BytesIO() as buffer_2:
            fits.writeto(buffer_2, image)
            with gzip.open(buffer_1, "wb") as gzipped_f:
                gzipped_f.write(buffer_2.getvalue())
        compressed = buffer_1.getvalue()

    return compressed
