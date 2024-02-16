"""
Module to calculate zero point by including a color term.
"""

import logging
from typing import Callable

import numpy as np
from astropy.table import Table
from scipy.odr import ODR, Model, RealData

from mirar.data import Image
from mirar.processors.photcal.zp_calculator.base_zp_calculator import (
    BaseZeroPointCalculator,
)

logger = logging.getLogger(__name__)


def line_func(theta: (float, float), x: np.ndarray) -> np.ndarray:
    """
    linear model to hand to scipy.odr
    Args:
        theta: slope and intercept of line to fit to data
        x: x data (color in our context)

    Returns:

    """
    slope, intercept = theta
    return (slope * x) + intercept


class ZPWithColorTermCalculator(
    BaseZeroPointCalculator
):  # pylint: disable=too-few-public-methods
    """
    Class to calculate zero point by including a color term. This models the data as

    ref_mag - img_mag = ZP + C * (ref_color).

    Note: this processor relies on scipy.odr,
    which requires all data points to have nonzero (x and y) error. This function
    removes those points, but it is up to the user to also consider removing these
    points from their reference catalogs to begin with (0 error in magnitude indicates
    unreliable photometry in PanStarrs, for example).

    Attributes:
        color_colnames_guess_generator: function that takes an image as input and
        returns three tuples. The first two tuples contain two strings each that are the
        column names of the reference catalog magnitudes and magnitude errors to use for
        the color term. The first string is the bluer band, the second is the redder
        band. The third tuple returns a first guess at the color and zero-point values.
    """

    def __init__(
        self,
        color_colnames_guess_generator: Callable[
            [Image], tuple[tuple[str, str], tuple[str, str], tuple[float, float]]
        ],
        reject_outliers: bool = True,
        num_stars_threshold: int = 5,
    ):
        self.color_colnames_guess_generator = color_colnames_guess_generator
        self.reject_outliers = reject_outliers
        self.num_stars_threshold = num_stars_threshold

    def calculate_zeropoint(  # pylint: disable=too-many-locals
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        (
            color_colnames,
            color_err_colnames,
            firstguess_color_zp,
        ) = self.color_colnames_guess_generator(image)
        colors = matched_ref_cat[color_colnames[0]] - matched_ref_cat[color_colnames[1]]

        for colname in colnames:
            y = matched_ref_cat["magnitude"] - matched_img_cat[colname]
            x = colors
            y_err = np.sqrt(
                matched_img_cat[colname.replace("MAG", "MAGERR")] ** 2
                + matched_ref_cat["magnitude_err"] ** 2
            )

            x_err = np.sqrt(
                matched_ref_cat[color_err_colnames[0]] ** 2
                + matched_ref_cat[color_err_colnames[1]] ** 2
            )

            # use scipy.odr to fit a line to data with x and y uncertainties
            # setup: remove sources with 0 uncertainty (or else scipy.odr won't work)
            zero_mask = (y_err == 0) | (x_err == 0)
            if np.sum(zero_mask) != 0:
                logger.debug(
                    f"Found {np.sum(zero_mask)} source(s) with zero reported "
                    f"uncertainty, removing them from calibrations."
                )
                x, y = x[~zero_mask], y[~zero_mask]
                x_err, y_err = x_err[~zero_mask], y_err[~zero_mask]

            line_model = Model(line_func)

            if self.reject_outliers:
                y_lo, y_up = np.percentile(y, [1, 99])
                outlier_mask = (y > y_lo) & (y < y_up)
                y, y_err, x, x_err = (
                    y[outlier_mask],
                    y_err[outlier_mask],
                    x[outlier_mask],
                    x_err[outlier_mask],
                )
                data = RealData(x, y, sx=x_err, sy=y_err)
                odr = ODR(data, line_model, beta0=firstguess_color_zp)
                out = odr.run()
                color, zero_point = out.beta
                y_pred = line_func((color, zero_point), x)
                y_residual = y - y_pred
                residual_rms = 0.5 * (
                    np.percentile(y_residual, 84.13) - np.percentile(y_residual, 15.86)
                )
                residual_outlier_mask = np.abs(y_residual) <= 4 * residual_rms

                if np.sum(residual_outlier_mask) < self.num_stars_threshold:
                    logger.warning(
                        f"Too few stars ({np.sum(residual_outlier_mask)}) "
                        f"to calculate zeropoint for {colname}. Res"
                    )
                y, y_err, x, x_err = (
                    y[residual_outlier_mask],
                    y_err[residual_outlier_mask],
                    x[residual_outlier_mask],
                    x_err[residual_outlier_mask],
                )

            if len(y) < self.num_stars_threshold:
                logger.warning(
                    f"Too few stars ({len(y)}) to calculate zeropoint for {colname}"
                )
                zero_point = -99.0
                zp_err = -99.0
                color = -99.0
                color_err = -99.0
            else:
                data = RealData(x, y, sx=x_err, sy=y_err)
                odr = ODR(data, line_model, beta0=firstguess_color_zp)
                out = odr.run()
                color, zero_point = out.beta
                color_err, zp_err = np.sqrt(np.diag(out.cov_beta))

            aperture = colname.split("_")[-1]
            image[f"ZP_{aperture}"] = zero_point
            image[f"ZP_{aperture}_std"] = zp_err
            image[f"ZP_{aperture}_nstars"] = len(y)
            image[f"C_{aperture}"] = color
            image[f"C_{aperture}_std"] = color_err

        return image
