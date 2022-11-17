import os

import astropy.io.fits
import numpy as np
from winterdrp.paths import base_raw_dir, raw_img_sub_dir
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import core_fields
import logging
from collections.abc import Callable
from winterdrp.io import open_fits
from glob import glob
from winterdrp.errors import ProcessorError, ImageNotFoundError
from winterdrp.data import Image, ImageBatch

logger = logging.getLogger(__name__)


class ImageLoader(BaseImageProcessor):
    """Processor to load raw images.
    """

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

    def __str__(self):
        return f"Processor to load images from the {self.input_sub_dir} subdirectory " \
               f"using the '{self.load_image.__name__}' function"

    def open_raw_image(
            self,
            path: str
    ) -> Image:

        data, header = self.load_image(path)

        for key in core_fields:
            if key not in header.keys():
                err = f"Essential key {key} not found in header. " \
                      f"Please add this field first. Available fields are: {list(header.keys())}"
                logger.error(err)
                raise KeyError(err)

        return Image(data.astype(np.float64), header)

    def open_raw_image_batch(
            self,
            paths: list
    ) -> ImageBatch:
        image_batch = ImageBatch()
        for path in paths:
            image = self.open_raw_image(path)
            image_batch.append(image)
        return image_batch

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:

        input_dir = os.path.join(
            self.input_img_dir,
            os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        return load_from_dir(input_dir, open_f=self.open_raw_image)


def load_from_dir(input_dir: str, open_f: Callable[[str], Image]) -> ImageBatch:

    img_list = sorted(glob(f'{input_dir}/*.fits'))

    logger.info(f"Loading from {input_dir}, with {len(img_list)} images")

    if len(img_list) < 1:
        err = f"No images found in {input_dir}. Please check path is correct!"
        logger.error(err)
        raise ImageNotFoundError(err)

    images = ImageBatch()

    for path in img_list:
        image = open_f(path)
        images.append(image)

    return images


# class RecentCalLoader(ImageLoader):
#
#     def _apply_to_images(
#             self,
#             images: list[np.ndarray],
#             headers: list[astropy.io.fits.Header],
#     ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
#
#         input_dir = os.path.join(
#             self.input_img_dir,
#             os.path.join(self.night_sub_dir, self.input_sub_dir)
#         )
#
#         return self.load_from_dir(input_dir)









    


