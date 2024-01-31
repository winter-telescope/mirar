"""
Module for subtracting the background from an image using Sextractor
"""

from pathlib import Path

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.io import open_fits
from mirar.paths import BASE_NAME_KEY
from mirar.processors import BaseImageProcessor
from mirar.processors.astromatic import Sextractor
from mirar.processors.astromatic.sextractor.sextractor import sextractor_checkimg_map
from mirar.processors.base_processor import PrerequisiteError, logger


class SextractorBkgSubtractor(BaseImageProcessor):
    """
    Processor to subtract the background from an image using Sextractor.
    This processor requires that Sextractor has been run previously with
    CHECKIMAGE_TYPE -BACKGROUND. This processor loads and returns the
    background-subtracted image.
    """

    base_key = "sextractorbkgsubtractor"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        new_batch = ImageBatch()
        for image in batch:
            if sextractor_checkimg_map["-BACKGROUND"] not in image.header:
                raise PrerequisiteError(
                    f"{sextractor_checkimg_map['-BACKGROUND']} key not found in image. "
                    f"Sextractor must be run with CHECKIMAGE_TYPE -BACKGROUND before "
                    "running this processor"
                )
            bkgsub_path = Path(image[sextractor_checkimg_map["-BACKGROUND"]])
            bkgsub_data, bkgsub_header = open_fits(bkgsub_path)

            # Mask the data with the original image's mask
            mask = np.isnan(image.get_data())
            bkgsub_data[mask] = np.nan

            # Update headers with any new keys that may have been added
            for key in image.header:
                if key not in bkgsub_header:
                    bkgsub_header[key] = image.header[key]
            # Copy over BASENAME
            bkgsub_header[BASE_NAME_KEY] = image.header[BASE_NAME_KEY]
            bkgsub_image = Image(data=bkgsub_data, header=bkgsub_header)
            new_batch.append(bkgsub_image)

            # Delete the original background file
            bkgsub_path.unlink()
        return new_batch

    def get_sextractor_module(self) -> Sextractor:
        """
        Get the Sextractor module from the preceding steps
        """
        mask = [isinstance(x, Sextractor) for x in self.preceding_steps]
        return np.array(self.preceding_steps)[mask][-1]

    def check_prerequisites(
        self,
    ):
        """
        Check that Sextractor has been run previously with CHECKIMAGE_TYPE -BACKGROUND
        """
        mask = [isinstance(x, Sextractor) for x in self.preceding_steps[-1:]]
        if np.sum(mask) < 1:
            err = (
                f"{self.__module__} requires running {Sextractor} immediately "
                f"preceding this processor as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise PrerequisiteError(err)

        sextractor_checkimg_types = self.get_sextractor_module().checkimage_type
        if "-BACKGROUND" not in sextractor_checkimg_types:
            err = (
                f"{self.__module__} requires that Sextractor be run with "
                f"CHECKIMAGE_TYPE -BACKGROUND. However, the following "
                f"CHECKIMAGE_TYPEs were found: {sextractor_checkimg_types}."
            )
            logger.error(err)
            raise PrerequisiteError(err)
