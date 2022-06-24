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


def split_images_into_batches(
        images: list[np.ndarray],
        headers: list[astropy.io.fits.Header],
        split_key: str | list[str]
) -> list[list[list[np.ndarray], list[astropy.io.fits.Header]]]:

    if isinstance(split_key, str):
        split_key = [split_key]

    groups = dict()

    for i, header in enumerate(headers):
        uid = []

        for key in split_key:
            uid.append(str(header[key]))

        uid = "_".join(uid)

        if uid not in groups.keys():
            groups[uid] = [[images[i]], [header]]
        else:
            groups[uid][0] += [images[i]]
            groups[uid][1] += [header]

    return [x for x in groups.values()]


class ImageBatcher(BaseProcessor):

    base_key = "batch"

    def __init__(
            self,
            split_key: str | list[str],
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.split_key = split_key

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        return images, headers

    def update_batches(
        self,
        batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]]
    ) -> list[list[list[np.ndarray], list[astropy.io.fits.header]]]:

        new_batches = []

        for [images, headers] in batches:
            new_batches += split_images_into_batches(images, headers, split_key=self.split_key)

        return new_batches

