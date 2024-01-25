"""
Module for calculating zero points
"""
import logging
import warnings
from typing import Callable

import numpy as np
from astropy.io.fits.verify import VerifyWarning
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.table import Table
from scipy.optimize import curve_fit

from mirar.data import Image
from mirar.processors.photcal import PhotometryCalculationError

logger = logging.getLogger(__name__)


class BaseZeroPointCalculator:
    """
    Base class for actually calculating zero point calculators
    """

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        """
        Function to calculate zero point from two catalogs and add the zeropoint
        information to the image header
        Args:
            matched_ref_cat: Reference catalog table
            matched_img_cat: Catalog of sources
            image: Image object
            colnames: List of column names from the image catalog to use for
            calculating zero point. The reference catalog is assumed to have the
            "magnitude" column, as it comes from mirar.catalog.base_catalog.BaseCatalog
        Returns:
            Image object with zero point information added to the header.
        """
        raise NotImplementedError


class OutlierRejectionZPCalculator(BaseZeroPointCalculator):
    """
    Class to calculate zero point using outlier rejection
    Attributes:
        num_stars_threshold: int to use as minimum number of stars
        outlier_rejection_threshold: float or list of floats to use as number of sigmas
    """

    def __init__(
        self,
        num_stars_threshold: int = 5,
        outlier_rejection_threshold: float | list[float] = 3.0,
    ):
        self.num_stars_threshold = num_stars_threshold
        self.outlier_rejection_threshold = outlier_rejection_threshold
        if isinstance(outlier_rejection_threshold, float):
            self.outlier_rejection_threshold = [outlier_rejection_threshold]
        self.outlier_rejection_threshold = np.sort(self.outlier_rejection_threshold)

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        zeropoints = []

        for colname in colnames:
            offsets = np.ma.array(
                matched_ref_cat["magnitude"] - matched_img_cat[colname]
            )
            for outlier_thresh in self.outlier_rejection_threshold:
                cl_offset = sigma_clip(offsets, sigma=outlier_thresh)
                num_stars = np.sum(np.invert(cl_offset.mask))

                zp_mean, zp_med, zp_std = sigma_clipped_stats(
                    offsets, sigma=outlier_thresh
                )

                if num_stars > self.num_stars_threshold:
                    break

            check = [np.isnan(x) for x in [zp_mean, zp_med, zp_std]]
            if np.sum(check) > 0:
                err = (
                    f"Error with nan when calculating sigma stats: \n "
                    f"mean: {zp_mean}, median: {zp_med}, std: {zp_std}"
                )
                logger.error(err)
                raise PhotometryCalculationError(err)

            aperture = colname.split("_")[-1]
            zero_dict = {
                "diameter": aperture,
                "zp_mean": zp_mean,
                "zp_median": zp_med,
                "zp_std": zp_std,
                "nstars": num_stars,
                "mag_cat": matched_ref_cat["magnitude"][np.invert(cl_offset.mask)],
                "mag_apers": matched_img_cat[colname][np.invert(cl_offset.mask)],
            }
            zeropoints.append(zero_dict)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore", category=VerifyWarning)

            for zpvals in zeropoints:
                image[f"ZP_{zpvals['diameter']}"] = zpvals["zp_mean"]
                image[f"ZP_{zpvals['diameter']}_std"] = zpvals["zp_std"]
                image[f"ZP_{zpvals['diameter']}_nstars"] = zpvals["nstars"]

        return image


def line_model(data, slope, intercept):
    """linear model to hand to scipy curve_fit"""
    return slope * data + intercept


class ZPWithColorTermCalculator(BaseZeroPointCalculator):
    """
    Class to calculate zero point by including a color term. This models the data as

    ref_mag - img_mag = ZP + C * (ref_color)

    Attributes:
        color_colnames_generator: function that takes an image as input and returns
        a list containing two strings that are the column names of the reference catalog
        to use for the color term. The first string is the bluer band, the second is the
         redder band.
    """

    def __init__(
        self,
        color_colnames_generator: Callable[[Image], [str, str]],
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
            y = matched_img_cat[colname] - matched_ref_cat["magnitude"]
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
