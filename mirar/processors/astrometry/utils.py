"""
Module containing processors to add astrometry headers to images
"""
import logging
import os

from astropy.io import fits

from mirar.data import ImageBatch
from mirar.paths import ASTROMETRY_FILE_KEY, BASE_NAME_KEY, get_astrometry_keys
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class AstrometryFromFile(BaseImageProcessor):
    """
    Processor to add astrometry headers to images from file.
    """

    base_key = "astrometry_from_file"

    def __init__(self, astrometry_file_key: str = ASTROMETRY_FILE_KEY):
        super().__init__()
        self.astrometry_file_key = astrometry_file_key

    def __str__(self) -> str:
        return "Processor to add astrometry headers to images from file."

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        astrometry_keys = get_astrometry_keys()
        new_batch = ImageBatch()
        for image in batch:
            astrometry_file = image[self.astrometry_file_key]

            if not os.path.exists(astrometry_file):
                raise FileNotFoundError(
                    f"Could not find astrometry file "
                    f"{astrometry_file}. Are you for sure running "
                    f"scamp with cache=True?"
                )
            for key in astrometry_keys:
                if key in image.header.keys():
                    logger.debug(f"Removing {key} from {image[BASE_NAME_KEY]}")
                    del image.header[key]

            logger.info(
                f"Adding astrometry headers from {astrometry_file} "
                f"to {image[BASE_NAME_KEY]}"
            )

            with open(astrometry_file, "r") as f:
                header_data = f.read()
                # Scamp v 2.10.0 writes this annoying character in a comment
                header_data = header_data.replace("é", "e")

            astrometry_header = fits.Header.fromstring(header_data, sep="\n")
            for k in astrometry_header:
                if (k == "HISTORY") | (k == "COMMENT"):
                    continue

                logger.debug(f"Adding {k} to {astrometry_header[k]}")
                image.header.append(
                    (k, astrometry_header[k], astrometry_header.comments[k])
                )
            new_batch.append(image)
        return new_batch
