"""
Core ZOGY algorithm implementation in Python.
############################################################
# Python implementation of ZOGY image subtraction algorithm
# See Zackay, Ofek, and Gal-Yam 2016 for details
# http://arxiv.org/abs/1601.02655
# SBC - 6 July 2016
# FJM - 20 October 2016
# SBC - 28 July 2017
# RDS - 30 October 2022
############################################################
"""

import logging
from pathlib import Path

import numpy as np
import pyfftw
import pyfftw.interfaces.numpy_fft as fft
from astropy.io import fits

logger = logging.getLogger(__name__)

pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(1.0)


def pyzogy(
    new_data: np.ndarray,
    ref_data: np.ndarray,
    new_psf_path: str | Path,
    ref_psf_path: str | Path,
    new_sigma_path: str | Path,
    ref_sigma_path: str | Path,
    new_avg_unc: float,
    ref_avg_unc: float,
    dx: float = 0.25,
    dy: float = 0.25,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Python implementation of ZOGY image subtraction algorithm.
    As per Frank's instructions, will assume images have been aligned,
    background subtracted, and gain-matched.

    Arguments:
    new_image_path: New image (filename)
    ref_image_path: Reference image (filename)
    new_psf_path: PSF of New image (filename)
    ref_psf_path: PSF or Reference image (filename)
    new_sigma_path: 2D Uncertainty (sigma) of New image (filename)
    ref_sigma_path: 2D Uncertainty (sigma) of Reference image (filename)
    new_avg_unc: Average uncertainty (sigma) of New image
    ref_avg_unc: Average uncertainty (sigma) of Reference image
    dx: Astrometric uncertainty (sigma) in x coordinate
    dy: Astrometric uncertainty (sigma) in y coordinate

    Returns:
    diff: Subtracted image
    diff_psf: PSF of subtracted image
    s_corr: Corrected subtracted image
    """

    # Set nans to zero in new and ref images
    new_nanmask = np.isnan(new_data)
    ref_nanmask = np.isnan(ref_data)

    new_data[new_nanmask] = 0.0
    ref_data[ref_nanmask] = 0.0

    # Load the PSFs into memory
    with fits.open(new_psf_path) as img_psf_f:
        new_psf = img_psf_f[0].data  # pylint: disable=no-member

    with fits.open(ref_psf_path) as ref_psf_f:
        ref_psf = ref_psf_f[0].data  # pylint: disable=no-member

    logger.debug(
        f"Max of small PSF is "
        f"{np.unravel_index(np.argmax(new_psf, axis=None), new_psf.shape)}"
    )

    # Place PSF at center of image with same size as new / reference
    new_psf_big = np.zeros(new_data.shape)
    ref_psf_big = np.zeros(ref_data.shape)

    y_min = new_data.shape[0] // 2 - new_psf.shape[0] // 2
    y_max = new_data.shape[0] // 2 + new_psf.shape[0] // 2 + 1
    x_min = new_data.shape[1] // 2 - new_psf.shape[1] // 2
    x_max = new_data.shape[1] // 2 + new_psf.shape[1] // 2 + 1

    new_psf_big[y_min:y_max, x_min:x_max] = new_psf
    ref_psf_big[y_min:y_max, x_min:x_max] = ref_psf

    logger.debug(
        f"Max of big PSF is "
        f"{np.unravel_index(np.argmax(new_psf_big, axis=None), new_psf_big.shape)}"
    )

    # Shift the PSF to the origin, so that it will not introduce a shift
    new_psf_big = fft.fftshift(new_psf_big)
    ref_psf_big = fft.fftshift(ref_psf_big)

    logger.debug(
        f"Max of big PSF shift is "
        f"{np.unravel_index(np.argmax(new_psf_big, axis=None), new_psf_big.shape)}"
    )

    # Take all the Fourier Transforms
    new_hat = fft.fft2(new_data)
    ref_hat = fft.fft2(ref_data)

    new_psf_hat = fft.fft2(new_psf_big)
    ref_psf_hat = fft.fft2(ref_psf_big)

    # Fourier Transform of Difference Image (Equation 13)
    diff_hat_numerator = ref_psf_hat * new_hat - new_psf_hat * ref_hat
    diff_hat_denominator = np.sqrt(
        new_avg_unc**2 * np.abs(ref_psf_hat**2)
        + ref_avg_unc**2 * np.abs(new_psf_hat**2)
        + 1e-8
    )

    diff_hat = diff_hat_numerator / diff_hat_denominator
    # Flux-based zero point (Equation 15)
    flux_zero_point = 1.0 / np.sqrt(new_avg_unc**2 + ref_avg_unc**2)

    # Difference Image
    diff = np.real(fft.ifft2(diff_hat)) / flux_zero_point
    # Fourier Transform of PSF of Subtraction Image (Equation 14)
    diff_hat_psf = ref_psf_hat * new_psf_hat / flux_zero_point / diff_hat_denominator

    # PSF of Subtraction Image
    diff_psf = np.real(fft.ifft2(diff_hat_psf))
    diff_psf = fft.ifftshift(diff_psf)
    diff_psf = diff_psf[y_min:y_max, x_min:x_max]
    logger.debug(
        f"Max of diff PSF is "
        f"{np.unravel_index(np.argmax(diff_psf, axis=None), diff_psf.shape)}"
    )

    # Fourier Transform of Score Image (Equation 17)
    score_hat = flux_zero_point * diff_hat * np.conj(diff_hat_psf)

    # Score Image
    score = np.real(fft.ifft2(score_hat))

    # Now start calculating Scorr matrix (including all noise terms)

    # Start out with source noise
    # Load the sigma images into memory
    with fits.open(new_sigma_path) as img_sigma_f:
        new_sigma = img_sigma_f[0].data  # pylint: disable=no-member

    with fits.open(ref_sigma_path) as ref_sigma_f:
        ref_sigma = ref_sigma_f[0].data  # pylint: disable=no-member

    new_sigma[new_nanmask] = 0.0
    ref_sigma[ref_nanmask] = 0.0

    # Sigma to variance
    new_variance = new_sigma**2
    ref_variance = ref_sigma**2

    # Fourier Transform of variance images
    new_variance_hat = fft.fft2(new_variance)
    ref_variance_hat = fft.fft2(ref_variance)

    # Equation 28
    k_r_hat = (
        np.conj(ref_psf_hat) * np.abs(new_psf_hat**2) / (diff_hat_denominator**2)
    )
    k_r = np.real(fft.ifft2(k_r_hat))

    # Equation 29
    k_n_hat = (
        np.conj(new_psf_hat) * np.abs(ref_psf_hat**2) / (diff_hat_denominator**2)
    )
    k_n = np.real(fft.ifft2(k_n_hat))

    # Noise in New Image: Equation 26
    new_noise = np.real(fft.ifft2(new_variance_hat * fft.fft2(k_n**2)))

    # Noise in Reference Image: Equation 27
    ref_noise = np.real(fft.ifft2(ref_variance_hat * fft.fft2(k_r**2)))

    # Astrometric Noise
    # Equation 31
    new_sigma = np.real(fft.ifft2(k_n_hat * new_hat))
    dsn_dx = new_sigma - np.roll(new_sigma, 1, axis=1)
    dsn_dy = new_sigma - np.roll(new_sigma, 1, axis=0)

    # Equation 30
    v_ast_s_n = dx**2 * dsn_dx**2 + dy**2 * dsn_dy**2

    # Equation 33
    ref_sigma = np.real(fft.ifft2(k_r_hat * ref_hat))
    dsr_dx = ref_sigma - np.roll(ref_sigma, 1, axis=1)
    dsr_dy = ref_sigma - np.roll(ref_sigma, 1, axis=0)

    # Equation 32
    v_ast_s_r = dx**2 * dsr_dx**2 + dy**2 * dsr_dy**2

    # Calculate Scorr
    s_corr = score / np.sqrt(new_noise + ref_noise + v_ast_s_n + v_ast_s_r)

    # Set back nans before returning
    diff[new_nanmask | ref_nanmask] = np.nan
    s_corr[new_nanmask | ref_nanmask] = np.nan

    return diff, diff_psf, s_corr
