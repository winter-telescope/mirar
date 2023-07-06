"""
Module for getting the component images used to make a swarp stack
"""
import logging
from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from mirar.data import Image, ImageBatch
from mirar.io import open_fits
from mirar.paths import STACKED_COMPONENT_IMAGES_KEY
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.base_processor import BaseImageProcessor
from mirar.processors.utils.image_saver import ImageSaver

logger = logging.getLogger(__name__)


class ReloadSwarpComponentImages(BaseImageProcessor):
    """
    Get the component images used to make a swarp stack
    """

    base_key = "swarp_component_images"

    def __init__(
        self,
        load_image: Callable[[str], [np.ndarray, Header]] = open_fits,
        header_key=STACKED_COMPONENT_IMAGES_KEY,
        copy_header_keys: str | list[str] = None,
    ):
        super().__init__()
        self.load_image = load_image
        self.header_key = header_key

        if isinstance(copy_header_keys, str):
            copy_header_keys = [copy_header_keys]
        self.copy_header_keys = copy_header_keys

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        if len(batch) > 1:
            raise NotImplementedError(
                "GetSwarpComponentImages only works on a batch containing a "
                "single images. Consider adding an ImageDebatcher before "
                "this processor."
            )
        component_batch = ImageBatch()
        image = batch[0]
        component_images_list = image[self.header_key].split(",")

        for component_image_path in component_images_list:
            if not Path(component_image_path).exists():
                raise FileNotFoundError(
                    f"Component image {component_image_path} not found. "
                    f"Are you sure it was saved using ImageSaver to this path just "
                    f"before the Swarp processor that stacked it?"
                )
            component_data, component_header = self.load_image(component_image_path)
            component_image = Image(component_data, component_header)
            if self.copy_header_keys is not None:
                for key in self.copy_header_keys:
                    if key in image.keys():
                        component_image[key] = image[key]
            component_batch.append(component_image)
        logger.info(f"Loaded {len(component_batch)} component images")
        return component_batch

    def check_prerequisites(
        self,
    ):
        mask = np.array([isinstance(x, Swarp) for x in self.preceding_steps])
        if np.sum(mask) == 0:
            err = (
                f"{self.__module__} requires {Swarp} as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise ValueError(err)

        index = np.argmax(mask)

        preceding_step = self.preceding_steps[index - 1]

        if not isinstance(preceding_step, ImageSaver):
            err = (
                f"{self.__module__} requires an {ImageSaver} to be used to save the "
                f"component images immediately before {Swarp} is run. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise ValueError(err)
