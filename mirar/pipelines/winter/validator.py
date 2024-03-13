"""
Module to reject images based on quality criteria
"""

import logging

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors.exceptions import ProcessorError
from mirar.paths import EXPTIME_KEY
from mirar.processors.astrometry.validate import PoorAstrometryError, PoorFWHMError
from mirar.processors.split import SUB_ID_KEY

logger = logging.getLogger(__name__)


class TooManyMaskedPixelsError(ProcessorError):
    """
    Error for when too many pixels in an image are masked
    """


class CondensationError(ProcessorError):
    """
    Error for when an image is affected by condensed
    """


class DarkOverSubtractionError(ProcessorError):
    """
    Error for when an image is affected by dark over-subtraction
    """


def masked_images_rejector(batch: ImageBatch) -> ImageBatch:
    """
    Rejects images with too many masked pixels
    """
    assert len(batch) == 1

    subdet_nan_limits = {
        1: 0.6,
        2: 0.6,
        3: 0.6,
        4: 0.6,
        5: 0.6,
        6: 0.6,
    }

    for image in batch:
        mask = image.get_mask()
        frac_masked = np.sum(~mask) / mask.size
        if frac_masked > subdet_nan_limits[image.header[SUB_ID_KEY]]:
            raise TooManyMaskedPixelsError(
                f"Fraction of masked pixels ({frac_masked}) is above threshold "
                f"{subdet_nan_limits[image.header[SUB_ID_KEY]]}"
            )
    return batch


def poor_astrometric_quality_rejector(batch: ImageBatch) -> ImageBatch:
    """
    Rejects images with poor astrometric quality
    1. Rejects images with SCAMP astrometric-reference RMS error above 0.3 arcsec
    2. Rejects images with median astrometric uncertainty (by comparing to Gaia2MASS)
    above 1.0 arcsec
    3. Rejects images with median FWHM above 6.0 arcsec
    """
    astrometric_unc_threshold_arcsec = 3.0
    fwhm_threshold_arcsec = 6.0
    for image in batch:
        if image["ASTUNC"] > astrometric_unc_threshold_arcsec / 3600:
            raise PoorAstrometryError(
                f"Uncertainty in astrometric solution from Scamp "
                f"({image['ASTUNC'] * 3600}) arcsec is above threshold "
                f"{astrometric_unc_threshold_arcsec} arcsec"
            )

        if image["FWHM_MED"] > fwhm_threshold_arcsec:
            raise PoorFWHMError(
                f"FWHM ({image['FWHM_MED']}) is above threshold"
                f" {fwhm_threshold_arcsec} arcsec."
            )
    return batch


def is_condensation_in_image(image: Image) -> bool:
    """
    Checks if a WINTER image is affected by condensation
    """
    data = image.get_data()
    header = image.get_header()
    vmedian = np.nanmedian(data, axis=1)
    x_inds = np.arange(len(vmedian))
    nanmask = np.invert(np.isnan(vmedian))
    wavmask = nanmask & (x_inds > 20) & (x_inds < 1070)
    boardid = header["BOARD_ID"]

    condensed = False
    if boardid in [1, 2, 3, 4, 5]:
        polydegs = np.polyfit(x=x_inds[wavmask], y=vmedian[wavmask], deg=1)
        interp1d_vals = np.polyval(polydegs, x_inds[wavmask])
        post_linear_trend_removal = vmedian[wavmask] / interp1d_vals
        pdeg2 = np.polyfit(x=x_inds[wavmask], y=post_linear_trend_removal, deg=2)
        interp2 = np.polyval(pdeg2, x_inds[wavmask])

        condensed = (pdeg2[0] > 0) & (np.ptp(interp2) > 0.1)

    return condensed


def winter_condensation_rejector(images: ImageBatch) -> ImageBatch:
    """
    Rejects images possibly affected by condensation
    """
    assert len(images) == 1
    for image in images:
        if is_condensation_in_image(image):
            raise CondensationError("Image is affected by condensation")
    return images


def winter_dark_oversubtraction_rejector(images: ImageBatch) -> ImageBatch:
    """
    Rejects images possibly affected by dark oversubtraction
    """
    assert len(images) == 1
    median_sky_counts_threshold_per_sec = 100.0 / 120.0
    for image in images:
        data = image.get_data()
        if np.nanmedian(data) < median_sky_counts_threshold_per_sec * image["EXPTIME"]:
            raise DarkOverSubtractionError(
                f"Dark-subtracted image has lower than expected median"
                f"counts for exposure time {image[EXPTIME_KEY]}."
                f"Threshold : {median_sky_counts_threshold_per_sec * image['EXPTIME']},"
                f" got: {np.nanmedian(data)}"
            )
    return images
