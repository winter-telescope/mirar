import astropy.io.fits
import numpy as np
import logging
from winterdrp.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


def filter_images(
        images: list[np.ndarray],
        headers: list[astropy.io.fits.Header],
        header_key: str = "target",
        reject_values: str | list[str] = "science",
) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

    # Enforce string in list for later matching
    if not isinstance(reject_values, list):
        reject_values = [str(reject_values)]
    else:
        reject_values = [str(x) for x in reject_values]

    passing_images = []
    passing_headers = []

    for i, header in enumerate(headers):
        if str(header[header_key]) not in reject_values:
            passing_images.append(images[i])
            passing_headers.append(header)

    return passing_images, passing_headers


class ImageRejector(BaseImageProcessor):

    base_key = "reject"

    def __init__(
            self,
            *args: tuple[str, str | list[str]],
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.rejects = args

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for (header_key, reject_values) in self.rejects:

            images, headers = filter_images(
                images,
                headers,
                header_key=header_key,
                reject_values=reject_values
            )

        return images, headers

    def update_batches(
        self,
        batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]]
    ) -> list[list[list[np.ndarray], list[astropy.io.fits.header]]]:
        # Remove empty batches
        return [x for x in batches if len(x[0]) > 0]



