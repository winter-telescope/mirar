"""
Module containing processors that mask cosmic rays
"""
import logging

import lacosmic
import numpy as np

from mirar.data import ImageBatch
from mirar.errors import NoncriticalProcessingError
from mirar.paths import EXPTIME_KEY
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class CRCleanError(NoncriticalProcessingError):
    """
    Error derived from running cosmic ray cleaning processor
    """


class LACosmicCleaner(BaseImageProcessor):
    """
    Processor to mask out cosmic rays.
    """

    base_key = "lacosmic"

    def __init__(
        self,
        contrast=2,
        cr_threshold=5,
        neighbor_threshold=0.3,
        error=None,
        background=None,
        effective_gain=None,
        readnoise=None,
        maxiter=4,
        border_mode="mirror",
        min_exptime=None,
        effective_gain_key=None,
        readnoise_key=None,
    ):
        super().__init__()
        self.contrast = contrast
        self.cr_threshold = cr_threshold
        self.neighbor_threshold = neighbor_threshold
        self.error = error
        self.background = background
        self.effective_gain = effective_gain
        self.readnoise = readnoise
        self.maxiter = maxiter
        self.border_mode = border_mode
        self.min_exptime = min_exptime
        self.effective_gain_key = effective_gain_key
        self.readnoise_key = readnoise_key

        if np.logical_and(self.effective_gain is None, self.effective_gain_key is None):
            err = (
                "LA Cosmic cleaner requires the gain. Please make sure you provide the "
                "gain value, or the header key to use for the gain. "
            )
            raise CRCleanError(err)

        if np.logical_and(self.readnoise is None, self.readnoise_key is None):
            err = (
                "LA Cosmic cleaner requires the read noise. Please make sure you "
                "provide the readnoise value, or the header key to use for the "
                "readnoise. "
            )
            raise CRCleanError(err)

        if np.logical_and(
            self.error is None,
            (
                (self.readnoise is None)
                & (self.effective_gain is None)
                & (self.effective_gain_key is None)
            ),
        ):
            err = (
                "LA Cosmic cleaner requires either a 2D error image, or the readnoise "
                "and the effective gain. Please make sure at least one is provided "
            )
            raise CRCleanError(err)

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        logger.info(f"Running cosmic ray cleaner on {len(batch)} images")

        for image in batch:
            run_crclean = True
            if self.min_exptime is not None:
                if image.header[EXPTIME_KEY] < self.min_exptime:
                    run_crclean = False
                    logger.warning(
                        "Exposure time is smaller than minimum specified, skipping "
                        "cosmic ray cleaner."
                    )

            if run_crclean:
                logger.debug("Running LACosmic")

                mask_data = ~image.get_mask().astype(bool)

                logger.debug(f"Mask data is : {mask_data}")
                effective_gain, readnoise = self.effective_gain, self.readnoise
                if effective_gain is None:
                    effective_gain = image.header[self.effective_gain_key]

                if readnoise is None:
                    readnoise = image.header[self.readnoise_key]
                logger.debug("Cleaning cosmic rays")
                cleaned_data, _ = lacosmic.lacosmic(
                    image.get_data(),
                    contrast=self.contrast,
                    cr_threshold=self.cr_threshold,
                    neighbor_threshold=self.neighbor_threshold,
                    error=self.error,
                    mask=mask_data,
                    background=self.background,
                    effective_gain=effective_gain,
                    readnoise=readnoise,
                    maxiter=self.maxiter,
                    border_mode=self.border_mode,
                )
                logger.debug("LACosmic finished cleaning")
                image.set_data(cleaned_data)

        return batch
