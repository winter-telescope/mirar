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


def line_func(theta, x):
    """linear model to hand to scipy.odr"""
    slope, intercept = theta
    return (slope * x) + intercept


class ZPWithColorTermCalculator(
    BaseZeroPointCalculator
):  # pylint: disable=too-few-public-methods
    """
    Class to calculate zero point by including a color term. This models the data as

    ref_mag - img_mag = ZP + C * (ref_color)

    Attributes:
        color_colnames_generator: function that takes an image as input and returns
        two lists containing two strings each that are the column names of the reference
        catalog magnitudes and magnitude errors to use for the color term. The first
        string is the bluer band, the second is the redder band.
    """

    def __init__(
        self,
        color_colnames_generator: Callable[
            [Image], list[list[str, str], list[str, str]]
        ],
    ):
        self.color_colnames_generator = color_colnames_generator

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
        ) = self.color_colnames_generator(image)
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
            ## setup: remove sources with 0 uncertainty (or else scipy.odr won't work)
            where_zero_y = np.where(np.array(y_err) == 0)[0]
            if len(where_zero_y) > 0:
                y = np.delete(y, where_zero_y)
                x = np.delete(x, where_zero_y)
                y_err = np.delete(y_err, where_zero_y)
                x_err = np.delete(x_err, where_zero_y)

            where_zero_x = np.where(np.array(x_err) == 0)[0]
            if len(where_zero_x) > 0:
                y = np.delete(y, where_zero_x)
                x = np.delete(x, where_zero_x)
                y_err = np.delete(y_err, where_zero_x)
                x_err = np.delete(x_err, where_zero_x)

            ## set up odr
            line_model = Model(line_func)
            data = RealData(x, y, sx=x_err, sy=y_err)
            odr = ODR(data, line_model, beta0=firstguess_color_zp)
            ## run the regression
            out = odr.run()
            color, zero_point = out.beta
            color_err, zp_err = np.sqrt(np.diag(out.cov_beta))

            aperture = colname.split("_")[-1]
            image[f"ZP_{aperture}"] = zero_point
            image[f"ZP_{aperture}_std"] = zp_err
            image[f"ZP_{aperture}_nstars"] = len(matched_ref_cat)
            image[f"C_{aperture}"] = color
            image[f"C_{aperture}_std"] = color_err

        return image
