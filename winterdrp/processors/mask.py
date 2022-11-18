import astropy.io.fits
import numpy as np
import logging
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.data import ImageBatch

logger = logging.getLogger(__name__)

# mask_value = -99.
mask_value = np.nan


class MaskPixels(BaseImageProcessor):

    base_key = "mask"

    def __init__(
            self,
            mask_path: str,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.mask = None
        self.mask_path = mask_path

    def __str__(self) -> str:
        return f"Processor to mask bad pixels using a pre-defined map: {self.mask_path}"

    def get_mask(self):
        if self.mask is None:
            self.mask = self.open_fits(self.mask_path)
        return self.mask

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:

        for i, image in enumerate(batch):
            data = image.get_data()
            mask = self.get_mask().get_data()
            mask = mask != 0
            data[mask] = mask_value
            image.set_data(data)

        return batch
