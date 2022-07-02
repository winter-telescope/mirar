import os

import astropy.io.fits
import numpy as np
from winterdrp.paths import base_raw_dir, raw_img_sub_dir
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.paths import core_fields
import logging
from collections.abc import Callable
from winterdrp.io import open_fits
from glob import glob

logger = logging.getLogger(__name__)


class ImageLoader(BaseProcessor):

    base_key = "load"

    def __init__(
            self,
            input_sub_dir: str = raw_img_sub_dir,
            input_img_dir: str = base_raw_dir,
            load_image: Callable = open_fits,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.input_sub_dir = input_sub_dir
        self.load_image = load_image
        self.input_img_dir = input_img_dir

    def open_raw_image(
            self,
            path: str
    ) -> tuple[np.array, astropy.io.fits.Header]:

        data, header = self.load_image(path)

        for key in core_fields:
            if key not in header.keys():
                err = f"Essential key {key} not found in header. " \
                      f"Please add this field first. Available fields are: {list(header.keys())}"
                logger.error(err)
                raise KeyError(err)

        return data.astype(np.float64), header

    def open_raw_image_batch(
            self,
            paths: list
    ) -> tuple[list, list]:

        images = []
        headers = []
        for path in paths:
            data, header = self.open_raw_image(path)
            images.append(data)
            headers.append(header)

        return images, headers

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        input_dir = os.path.join(
            self.input_img_dir,
            os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        img_list = glob(f'{input_dir}/*.fits')

        logger.info(f"Loading from {input_dir}, with {len(img_list)} images")

        new_images = []
        new_headers = []

        for path in img_list:
            img, header = self.open_raw_image(path)
            new_images.append(img)
            new_headers.append(header)

        return new_images, new_headers





    


