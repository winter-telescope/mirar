"""
Module to subtract a Fourier background model from a raw image.
"""
import numpy as np
from numpy import fft


def subtract_fourier_background_model(
    raw_data: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Subtract a Fourier background model from a raw image.
    """
    # USER-CONFIGURABLE PARAMETERS

    # Filter out stars from the original image before taking the FFT
    # (this suppresses broadband power in the FFT from point source visibilities
    # in the image space).  Middle parameter should always be 0.5, i.e. median.
    # Upper and lower can be configured but must be between 0-1.

    star_filter_params = [0.05, 0.5, 0.95]

    # This is used to filter out "point sources" in Fourier space corresponding to
    # high-power modes (i.e. electronic noise).  Setting the threshold too low could
    # filter out real sources, setting it too high will not filter the noise.  0.95
    # seemed like a sensible level in initial tests, must be between 0-1.
    fft_threshold = 0.95

    # there is some horizontal striping, take that out.
    vec = np.nanmedian(raw_data, axis=1)
    horizontal_stripes = np.outer(vec, np.ones(raw_data.shape[1]))

    data = raw_data - horizontal_stripes

    # Filter out stars (top 5%) and bad pixels (bottom 5%), set these pixels to the
    # median
    quantiles = np.quantile(data, star_filter_params)
    data[data > quantiles[2]] = quantiles[1]
    data[data < quantiles[0]] = quantiles[1]

    # Take the FFT, and find the brightest "point sources" in the transform.
    # These are discrete modes of electronic noise
    # This is a really brain-dead way of finding points, just take the brightest 5% of
    # modes.
    # Not clear how this will work on data with large extended sources.

    fourier_trans_data = fft.fft2(data)
    quantiles = np.quantile(np.abs(fourier_trans_data), [fft_threshold])

    # Make a noise-model using source-extractor, this is a bit of a hack that does not
    # work fully right now, but hopefully will soon.
    # run_sextractor('fourier.fits')
    # f_sex = fits.getdata('fourier.fits')
    # seg = fits.getdata('fourier.fits.seg.fits')
    # f_sex[seg == 0] = 0
    # model_sex = np.real(fft.ifft2(f_sex))

    # Make a low-noise model image that only includes the sharp modes
    fourier_trans_data[np.where(np.abs(fourier_trans_data) < quantiles[0])] = 0
    model = np.real(fft.ifft2(fourier_trans_data))

    # Note: filtered is the final output array, if you put this into a pipeline
    # this is the frame that you want to return.  It's possible that this should
    # have some checks to set more bad pixels to np.nan or make sure I haven't
    # mistakenly set some pixels incorrectly around the border.

    filtered_data = raw_data - horizontal_stripes - model
    # filtered[nans] = np.nan

    return filtered_data, model
