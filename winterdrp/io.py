from astropy.io import fits
import numpy as np
import astropy.io.fits


def create_fits(data, header):
    proc_hdu = fits.PrimaryHDU(data)
    if header is not None:
        proc_hdu.header = header  # Copy over the header from the raw file
    return proc_hdu


def save_to_path(
        data: np.ndarray,
        header: astropy.io.fits.Header,
        path: str,
        overwrite: bool = True
):
    img = create_fits(data, header=header)
    img.writeto(path, overwrite=overwrite)


def open_fits(
        path: str
) -> (np.array, astropy.io.fits.Header):
    with fits.open(path) as img:
        data = img[0].data
        header = img[0].header

    return data, header
