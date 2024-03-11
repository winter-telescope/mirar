"""
Module to calculate zero point by including a color term.
"""

import logging
from typing import Callable

import numpy as np
from astropy.table import MaskedColumn, Table
from scipy.odr import ODR, Model, RealData
from scipy.optimize import curve_fit

from mirar.data import Image
from mirar.processors.photcal.zp_calculator.base_zp_calculator import (
    BaseZeroPointCalculator,
)

logger = logging.getLogger(__name__)


def line_func_odr(theta: (float, float), x: np.ndarray) -> np.ndarray:
    """
    linear model to hand to scipy.odr
    Args:
        theta: slope and intercept of line to fit to data
        x: x data (color in our context)

    Returns:

    """
    slope, intercept = theta
    return (slope * x) + intercept


def line_func_curve_fit(data: np.ndarray, slope: float, intercept: float) -> np.ndarray:
    """linear model to hand to scipy curve_fit"""
    return slope * data + intercept


def solve_odr(
    x: np.ndarray,
    y: np.ndarray,
    x_err: np.ndarray,
    y_err: np.ndarray,
    firstguess_color_zp: tuple,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solve for the zero point and color term using scipy.odr.
    :param x: color
    :param y: ref_mag - img_mag
    :param x_err: uncertainty in color
    :param y_err: uncertainty in ref_mag - img_mag
    :param firstguess_color_zp: first guess at the color and zero-point values
    :return: best fit color and zero-point values, and their uncertainties
    """
    # setup: remove sources with 0 uncertainty (or else scipy.odr won't work)
    zero_mask = (y_err == 0) | (x_err == 0)
    if np.sum(zero_mask) != 0:
        logger.debug(
            f"Found {np.sum(zero_mask)} source(s) with zero reported "
            f"uncertainty, removing them from calibrations."
        )
        x, y = x[~zero_mask], y[~zero_mask]
        x_err, y_err = x_err[~zero_mask], y_err[~zero_mask]

    line_model = Model(line_func_odr)
    data = RealData(x, y, sx=x_err, sy=y_err)
    odr = ODR(data, line_model, beta0=firstguess_color_zp)
    out = odr.run()
    return out.beta, out.sd_beta


def solve_curve_fit(
    x: np.ndarray, y: np.ndarray, _, y_err: np.ndarray, firstguess_color_zp: tuple
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solve for the zero point and color term using scipy.curve_fit.
    :param x: color
    :param y: ref_mag - img_mag
    :param y_err: uncertainty in ref_mag - img_mag
    :param firstguess_color_zp: first guess at the color and zero-point values
    :return: best-fit color and zero-point values, and their uncertainties
    """
    popt, pcov = curve_fit(  # pylint: disable=unbalanced-tuple-unpacking
        f=line_func_curve_fit, xdata=x, ydata=y, sigma=y_err, p0=firstguess_color_zp
    )
    return popt, np.sqrt(np.diag(pcov))


class ZPWithColorTermCalculator(
    BaseZeroPointCalculator
):  # pylint: disable=too-few-public-methods
    """
    Class to calculate zero point by including a color term. This models the data as

    ref_mag - img_mag = ZP + C * (ref_color).

    Note: This processor can solve using either scipy.curve_fit or scipy.odr. scipy.odr
    is the default solver, and uses errors in both ref_color and ref_mag to fit the
    color term. scipy.curve_fit only uses errors in ref_mag to fit the color term. In
    both cases, the processor will remove any data points that are masked (missing
    entries) in the reference catalog. Additionally, scipy.odr requires all data points
    to have nonzero (x and y) error. The processor will removes those points, but it is
    up to the user to also consider removing these
    points from their reference catalogs to begin with (0 error in magnitude indicates
    unreliable photometry in PanStarrs, for example).

    Attributes:
        :param color_colnames_guess_generator: function that takes an image as input and
        returns three tuples. The first two tuples contain two strings each that are the
        column names of the reference catalog magnitudes and magnitude errors to use for
        the color term. The first string is the bluer band, the second is the redder
        band. The third tuple returns a first guess at the color and zero-point values.
        :param reject_outliers: whether to reject outliers when fitting the color term.
        :param num_stars_threshold: minimum number of stars to fit the color term.
        :param solver: solver to use for fitting the color term. Must be 'odr' or
        'curve_fit'. Defaults to 'odr'.
    """

    def __init__(
        self,
        color_colnames_guess_generator: Callable[
            [Image], tuple[tuple[str, str], tuple[str, str], tuple[float, float]]
        ],
        reject_outliers: bool = True,
        num_stars_threshold: int = 5,
        solver: str = "odr",
    ):
        self.color_colnames_guess_generator = color_colnames_guess_generator
        self.reject_outliers = reject_outliers
        self.num_stars_threshold = num_stars_threshold
        if solver not in ["odr", "curve_fit"]:
            raise ValueError("Solver must be 'odr' or 'curve_fit'")
        self.solver = solver
        if self.solver == "odr":
            self.solver_func = solve_odr
        else:
            self.solver_func = solve_curve_fit

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

        logger.debug(
            f"Start of calibrations: {len(matched_ref_cat)} cross-matched source(s). "
            f"Now for quality cuts..."
        )

        # If any of the columns in the reference catalog is a MaskedColumn, remove the
        # masked values
        for colname in color_colnames + color_err_colnames:
            if isinstance(matched_ref_cat[colname], MaskedColumn):
                mask = matched_ref_cat[colname].mask
                matched_ref_cat = matched_ref_cat[~mask]
                matched_img_cat = matched_img_cat[~mask]
                logger.debug(
                    f"Found {np.sum(mask)} source(s) with masked values in reference "
                    f"{colname} column, removing them from calibrations."
                )

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

            if self.reject_outliers:
                y_lo, y_up = np.percentile(y, [1, 99])
                outlier_mask = (y > y_lo) & (y < y_up)
                logger.debug(
                    f"Found {len(y) - np.sum(outlier_mask)} outlier source(s), "
                    f"removing them from calibrations."
                )
                y, y_err, x, x_err = (
                    y[outlier_mask],
                    y_err[outlier_mask],
                    x[outlier_mask],
                    x_err[outlier_mask],
                )

                (color, zero_point), _ = self.solver_func(
                    x, y, x_err, y_err, firstguess_color_zp
                )
                y_pred = line_func_odr((color, zero_point), x)
                y_residual = y - y_pred
                residual_rms = 0.5 * (
                    np.percentile(y_residual, 84.13) - np.percentile(y_residual, 15.86)
                )
                residual_outlier_mask = np.abs(y_residual) <= 4 * residual_rms

                if np.sum(residual_outlier_mask) < self.num_stars_threshold:
                    logger.warning(
                        f"Too few stars ({np.sum(residual_outlier_mask)}) "
                        f"to calculate zeropoint for {colname}."
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
                logger.debug(
                    f"End of calibrations: {len(y)} sources pass quality cuts."
                )
                (color, zero_point), (color_err, zp_err) = self.solver_func(
                    x, y, x_err, y_err, firstguess_color_zp
                )

            aperture = colname.split("_")[-1]
            image[f"ZP_{aperture}"] = zero_point
            image[f"ZP_{aperture}_std"] = zp_err
            image[f"ZP_{aperture}_nstars"] = len(y)
            image[f"C_{aperture}"] = color
            image[f"C_{aperture}_std"] = color_err

        return image
