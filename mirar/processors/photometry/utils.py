"""
Module with utis for photometry
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib.patches import Circle
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry

from mirar.data import Image
from mirar.errors import ProcessorError
from mirar.paths import GAIN_KEY

logger = logging.getLogger(__name__)


class CutoutError(ProcessorError):
    """
    Error raised when cutout generation fails
    """


def make_cutouts(
    image_paths: Path | list[Path], position: tuple, half_size: int
) -> list[np.array]:
    """
    Function to make cutouts
    Args:
        :param: image_paths: Path or list of paths to the images
        :param: position: (x,y) coordinates of the center of the cutouts
        :param: half_size: half_size of the square cutouts

    Returns:
        :return: cutout_list: list of 2D numpy arrays
    """
    if not isinstance(image_paths, list):
        image_paths = [image_paths]

    cutout_list = []
    for image_path in image_paths:
        data = fits.getdata(image_path)
        y_image_size, x_image_size = np.shape(data)
        x, y = position
        logger.debug(f"Cutout parameters {x},{y},{np.shape(data)}")

        if x < 0 or x > x_image_size or y < 0 or y > y_image_size:
            raise CutoutError(
                f"Cutout position {x},{y} is outside the image {image_path}"
            )

        if x < half_size:
            xmin, xmax = 0, x + half_size + 1
            x_pad_width = (half_size - x, 0)

        elif x + half_size + 1 > x_image_size:
            xmin, xmax = x - half_size, x_image_size
            x_pad_width = (0, (half_size + x + 1) - x_image_size)
        else:
            xmin, xmax = x - half_size, x + half_size + 1
            x_pad_width = (0, 0)

        if y < half_size:
            ymin, ymax = 0, y + half_size + 1
            y_pad_width = (half_size - y, 0)

        elif y + half_size + 1 > y_image_size:
            ymin, ymax = y - half_size, y_image_size
            y_pad_width = (0, (half_size + y + 1) - y_image_size)
        else:
            ymin, ymax = y - half_size, y + half_size + 1
            y_pad_width = (0, 0)

        cutout = data[ymin:ymax, xmin:xmax]
        cutout = np.pad(cutout, (y_pad_width, x_pad_width), "constant")
        cutout_list.append(cutout)

        assert cutout.shape == (2 * half_size + 1, 2 * half_size + 1), (
            f"Cutout shape is not correct for x={x} and y={y} and "
            f"half_size={half_size}. Shape is {cutout.shape} Data shape is {data.shape}"
        )
    return cutout_list


def psf_photometry(
    image_cutout: np.array,
    image_unc_cutout: np.array,
    psfmodels: np.array,
    psf_filename: str,
) -> tuple[float, float, float, float, float]:
    """
    Function to perform PSF photometry
    Args:
        :param: image_cutout: 2D numpy array of the image cutout
        :param: image_unc_cutout: 2D numpy array of the image uncertainty cutout
        :param: psfmodels: 3D numpy array of the PSF models

    Returns:
        :return: psf_fluxes: PSF flux
        :return: psf_flux_uncs: PSF flux uncertainty
        :return: chi2s: chi2 value
        :return xshifts: xshift required to match PSF to the source
        :return yshifts: yshift required to match PSF to the source
    """
    numpsfmodels = psfmodels.shape[2]

    chi2s, psf_fluxes, psf_flux_uncs = [], [], []
    for ind in range(numpsfmodels):
        psfmodel = psfmodels[:, :, ind]
        psf_flux = np.nansum(psfmodel * image_cutout) / np.nansum(np.square(psfmodel))
        psf_flux_unc = np.sqrt(
            np.nansum(np.square(psfmodel) * np.square(image_unc_cutout))
        ) / np.nansum(np.square(psfmodel))
        deg_freedom = np.size(image_cutout) - 1
        chi2 = (
            np.nansum(
                np.square(image_cutout - psfmodel * psf_flux)
                / np.square(image_unc_cutout)
            )
            / deg_freedom
        )

        psf_fluxes.append(psf_flux)
        psf_flux_uncs.append(psf_flux_unc)
        chi2s.append(chi2)

    minchi2_ind = np.argmin(chi2s)
    minchi2 = np.min(chi2s)
    best_fit_psf_flux = psf_fluxes[minchi2_ind]
    best_fit_psf_fluxunc = psf_flux_uncs[minchi2_ind]

    best_fit_psfmodel = psfmodels[:, :, minchi2_ind]
    ys_cen, xs_cen = np.where(best_fit_psfmodel == np.max(best_fit_psfmodel))

    psf = fits.getdata(psf_filename)
    y_cen, x_cen = np.where(psf == np.max(psf))
    yshift = ys_cen[0] - y_cen
    xshift = xs_cen[0] - x_cen

    return best_fit_psf_flux, best_fit_psf_fluxunc, minchi2, xshift, yshift


def make_psf_shifted_array(
    psf_filename: str, cutout_size_psf_phot: int = 20, pad_psf_size: int = 60
):
    """
    Function to make a shifted array from a PSF model
    Args:
        psf_filename:
        cutout_size_psf_phot:
        pad_psf_size:

    Returns:

    """
    psf = fits.getdata(psf_filename)
    normpsf = psf / np.sum(psf)
    ngrid = 81
    x_crds = np.linspace(-4, 4, 9)
    grid_x, grid_y = np.meshgrid(x_crds, x_crds)
    grid_x = np.ndarray.flatten(grid_x)
    grid_y = np.ndarray.flatten(grid_y)

    psfmodel_half_size = int(psf.shape[0] / 2)
    padpsfs = np.zeros((pad_psf_size, pad_psf_size, ngrid))
    pad_loind = int((pad_psf_size - 2 * psfmodel_half_size) / 2)
    pad_upind = pad_psf_size - pad_loind + 1
    for i in range(ngrid):
        padpsfs[
            int(pad_loind + grid_y[i]) : int(pad_upind + grid_y[i]),
            int(pad_loind + grid_x[i]) : int(pad_upind + grid_x[i]),
            i,
        ] = normpsf

    unshifted_ind = int(ngrid / 2)
    normpsfmax = np.max(normpsf)
    xcen_1, xcen_2 = np.where(padpsfs[:, :, unshifted_ind] == normpsfmax)
    xcen_1 = int(xcen_1)
    xcen_2 = int(xcen_2)

    psfmodels = padpsfs[
        xcen_1 - cutout_size_psf_phot : xcen_1 + cutout_size_psf_phot + 1,
        xcen_2 - cutout_size_psf_phot : xcen_2 + cutout_size_psf_phot + 1,
    ]
    return psfmodels


def aper_photometry(
    image_cutout: np.ndarray,
    image_unc_cutout: np.ndarray,
    aper_diameter: float,
    bkg_in_diameter: float,
    bkg_out_diameter: float,
    plot: bool = False,
    plotfilename: str = None,
):
    """
    Perform aperture photometry
    Args:
        :param: image_cutout: np.ndarray of the image cutout
        :param: image_unc_cutout: np.ndarray of the image uncertainty cutout
        :param: aper_diameter: aperture diameter in pixels
        :param: bkg_in_diameter: inner background annulus diameter in pixels
        :param: bkg_out_diameter: outer background annulus diameter in pixels
        :param: plot: whether to plot the cutout
        :param: plotfilename: filename to save plot to

    Returns:
        :return: aperture flux, aperture flux uncertainty
    """
    x_crd, y_crd = int(image_cutout.shape[0] / 2), int(image_cutout.shape[1] / 2)
    image_cutout_mask = np.isnan(image_cutout)
    if plot:
        if plotfilename is None:
            raise ValueError("Please provide a filename to save the plot to.")
        fig, plot_ax = plt.subplots()
        mean, std = np.nanmean(image_cutout), np.nanstd(image_cutout)
        plot_ax.imshow(
            image_cutout,
            interpolation="nearest",
            cmap="gray",
            vmin=mean - std,
            vmax=mean + 10 * std,
            origin="lower",
        )
        # c = Circle(xy=(x_img, y_img),radius=15)

        circle = Circle(xy=(x_crd, y_crd), radius=aper_diameter / 2)
        bkg_inner_annulus = Circle(xy=(x_crd, y_crd), radius=bkg_in_diameter / 2)
        bkg_outer_annulus = Circle(xy=(x_crd, y_crd), radius=bkg_out_diameter / 2)
        circle.set_facecolor("none")
        circle.set_edgecolor("red")
        bkg_inner_annulus.set_facecolor("none")
        bkg_inner_annulus.set_edgecolor("red")
        bkg_outer_annulus.set_facecolor("none")
        bkg_outer_annulus.set_edgecolor("red")
        plot_ax.add_artist(circle)
        plot_ax.add_artist(bkg_inner_annulus)
        plot_ax.add_artist(bkg_outer_annulus)
        plot_ax.set_xlim(x_crd - 30, x_crd + 30)
        plot_ax.set_ylim(y_crd - 30, y_crd + 30)
        plt.savefig(plotfilename)
        plt.close(fig)

    aperture = CircularAperture((x_crd, y_crd), r=aper_diameter)
    annulus_aperture = CircularAnnulus(
        (x_crd, y_crd), r_in=bkg_in_diameter / 2, r_out=bkg_out_diameter / 2
    )

    annulus_masks = annulus_aperture.to_mask(method="center")
    annulus_data = annulus_masks.multiply(image_cutout)
    mask = annulus_masks.data
    annulus_data_1d = annulus_data[mask > 0]
    _, bkg_median, _ = sigma_clipped_stats(annulus_data_1d, sigma=2, mask_value=np.nan)
    bkg = np.zeros(image_cutout.shape) + bkg_median
    # bkg_error = np.zeros(image_cutout.shape) + bkg_std

    aperture_mask = aperture.to_mask(method="center")
    aperture_unc_data = aperture_mask.multiply(image_unc_cutout)
    # effective_gain = header['GAIN']
    # error = calc_total_error(data, bkg_error, effective_gain)
    # phot_table = aperture_photometry(diff_cutout - bkg, aperture, error=error)
    # counts_err = phot_table['aperture_sum_err'][0]
    error = np.sqrt(np.nansum(aperture_unc_data**2))
    phot_table = aperture_photometry(
        data=image_cutout - bkg, apertures=aperture, mask=image_cutout_mask
    )
    counts = phot_table["aperture_sum"][0]
    counts_err = error
    return counts, counts_err


def get_rms_image(image: Image) -> Image:
    """Get an RMS image from a regular image

    :param image: An :class:`~mirar.data.image_data.Image`
    :param rms: rms of the image
    :return: An RMS :class:`~mirar.data.image_data.Image`
    """
    image_data = image.get_data()
    image_data = image_data[np.invert(np.isnan(image_data))]
    rms = 0.5 * (
        np.percentile(image_data[image_data != 0.0], 84.13)
        - np.percentile(image_data[image_data != 0.0], 15.86)
    )
    gain = image[GAIN_KEY]
    poisson_noise = np.copy(image.get_data()) / gain
    poisson_noise[poisson_noise < 0] = 0
    rms_image = Image(data=np.sqrt(poisson_noise + rms**2), header=image.get_header())
    return rms_image


def get_mags_from_fluxes(
    flux_list: np.ndarray[float],
    fluxunc_list: np.ndarray[float],
    zeropoint: float,
    zeropoint_unc: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert fluxes and flux uncertainties to magnitudes and magnitude uncertainties

    :param flux_list: List of fluxes
    :param fluxunc_list: List of flux uncertainties
    :param zeropoint: Zeropoint
    :param zeropoint_unc: Zeropoint uncertainties
    :return: magnitudes, magnitude uncertainties
    """
    assert len(flux_list) == len(fluxunc_list)

    magnitudes = zeropoint - 2.5 * np.log10(flux_list)
    magnitudes_unc = 1.086 * fluxunc_list / flux_list

    magnitudes_unc = np.sqrt(magnitudes_unc**2 + zeropoint_unc**2)
    magnitudes = magnitudes.astype(object)
    magnitudes_unc = magnitudes_unc.astype(object)
    magnitudes[flux_list <= 0] = None
    magnitudes_unc[flux_list <= 0] = None
    return magnitudes, magnitudes_unc
