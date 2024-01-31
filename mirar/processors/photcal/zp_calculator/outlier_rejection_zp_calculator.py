"""
Module to calculate zero point using outlier rejection
"""

import logging
import warnings

import numpy as np
from astropy.io.fits.verify import VerifyWarning
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.table import Table

from mirar.data import Image
from mirar.processors.photcal.photcal_errors import PhotometryCalculationError
from mirar.processors.photcal.zp_calculator.base_zp_calculator import (
    BaseZeroPointCalculator,
)

logger = logging.getLogger(__name__)


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
