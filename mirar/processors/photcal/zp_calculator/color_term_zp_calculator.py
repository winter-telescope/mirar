"""
Module to calculate zero point by including a color term.
"""

import logging
from typing import Callable

import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

from mirar.data import Image
from mirar.processors.photcal.zp_calculator.base_zp_calculator import (
    BaseZeroPointCalculator,
)

logger = logging.getLogger(__name__)


def line_model(data, slope, intercept):
    """linear model to hand to scipy curve_fit"""
    return slope * data + intercept


class ZPWithColorTermCalculator(BaseZeroPointCalculator):
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

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        color_colnames = self.color_colnames_generator(image)
        colors = matched_ref_cat[color_colnames[0]] - matched_ref_cat[color_colnames[1]]

        for colname in colnames:
            y = matched_ref_cat["magnitude"] - matched_img_cat[colname]
            x = colors
            y_err = np.sqrt(
                matched_img_cat[colname.replace("MAG", "MAGERR")] ** 2
                + matched_ref_cat["magnitude_err"] ** 2
            )
            popt, pcov = curve_fit(f=line_model, xdata=x, ydata=y, sigma=y_err)
            color, zero_point = popt
            color_err, zp_err = np.sqrt(np.diag(pcov))

            aperture = colname.split("_")[-1]
            image[f"ZP_{aperture}"] = zero_point
            image[f"ZP_{aperture}_std"] = zp_err
            image[f"ZP_{aperture}_nstars"] = len(matched_ref_cat)
            image[f"C_{aperture}"] = color
            image[f"C_{aperture}_std"] = color_err

        return image
