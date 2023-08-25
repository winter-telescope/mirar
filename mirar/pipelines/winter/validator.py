"""
Module to reject images based on quality criteria
"""
import numpy as np

from mirar.data import ImageBatch
from mirar.errors.exceptions import ProcessorError
from mirar.processors.astrometry.validate import PoorAstrometryError, PoorFWHMError
from mirar.processors.split import SUB_ID_KEY


class TooManyMaskedPixelsError(ProcessorError):
    """
    Error for when too many pixels in an image are masked
    """


def masked_images_rejector(batch: ImageBatch) -> ImageBatch:
    """
    Rejects images with too many masked pixels
    """
    assert len(batch) == 1
    subdet_nan_limits = {
        1: 0.4,
        2: 0.5,
        3: 0.3,
        4: 0.3,
        5: 0.4,
        6: 0.5,
        7: 0.4,
        8: 0.3,
        9: 0.4,
        10: 0.3,
        11: 0.4,
        12: 0.4,
    }

    for image in batch:
        mask = image.get_mask()
        frac_masked = np.sum(~mask) / mask.size
        if frac_masked > subdet_nan_limits[image.header[SUB_ID_KEY]]:
            raise TooManyMaskedPixelsError(
                f"Fraction of masked pixels ({frac_masked}) is above threshold"
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
    scamp_astrrms_threshold_arcsec = 0.5
    astrometric_unc_threshold_arcsec = 1.0
    fwhm_threshold_arcsec = 6.0
    for image in batch:
        if (image["ASTRRMS1"] > scamp_astrrms_threshold_arcsec / 3600) | (
            image["ASTRRMS2"] > scamp_astrrms_threshold_arcsec / 3600
        ):
            raise PoorAstrometryError(
                f"RMS astrometric error from Scamp "
                f"({image['ASTRRMS1']}, {image['ASTRRMS2']}) is above threshold"
            )

        if image["ASTUNC"] > astrometric_unc_threshold_arcsec / 3600:
            raise PoorAstrometryError(
                f"Uncertainty in astrometric solution from Scamp "
                f"({image['ASTUNC']}) is above threshold"
            )

        if image["FWHM_MED"] > fwhm_threshold_arcsec:
            raise PoorFWHMError(f"FWHM ({image['FWHM_MED']}) is above threshold")
    return batch
