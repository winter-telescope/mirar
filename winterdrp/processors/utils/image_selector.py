import astropy.io.fits
import numpy as np
from winterdrp.processors.base_processor import BaseProcessor


def select_from_images(
        images: list[np.ndarray],
        headers: list[astropy.io.fits.Header],
        target: str | list[str]
) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

    if isinstance(target, str):
        target = [target]

    passing_images = []
    passing_headers = []

    for i, header in enumerate(headers):
        if header["target"] in target:
            passing_images.append(images[i])
            passing_headers.append(header)

    return passing_images, passing_headers


class ImageSelector(BaseProcessor):

    base_key = "select"

    def __init__(
            self,
            target: str | list[str] = "science",
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.target = target

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        return select_from_images(images, headers, target=self.target)

