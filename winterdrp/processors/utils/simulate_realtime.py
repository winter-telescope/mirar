import logging
import os
from pathlib import Path
from shutil import copy

import astropy.io.fits
import numpy as np

from winterdrp.data import ImageBatch
from winterdrp.paths import (
    LATEST_SAVE_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    base_output_dir,
    get_output_dir,
    get_output_path,
)
from winterdrp.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class RealtimeImageSimulator(BaseImageProcessor):

    base_key = "simrealtime"

    def __init__(
        self,
        input_img_dir: str,
        input_img_names: str | list[str],
        output_dir_name: str,
        output_dir: str = base_output_dir,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        if not isinstance(input_img_names, list):
            input_img_names = [input_img_names]

        self.input_img_names = input_img_names
        self.input_img_dir = Path(input_img_dir)
        self.output_dir_name = output_dir_name
        self.output_dir = output_dir

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:

        for image_name in self.input_img_names:

            img_path = self.input_img_dir.joinpath(image_name)

            output_path = Path(
                get_output_dir(
                    dir_root=self.output_dir_name,
                    sub_dir=self.night_sub_dir,
                    output_dir=self.output_dir,
                )
            )
            logger.debug(f"Copying {img_path} to {output_path}")

            if not output_path.exists():
                output_path.mkdir(parents=True)

            output_path = output_path.joinpath(img_path.name)

            if output_path.is_file():
                output_path.unlink()

            copy(img_path, output_path)

        return batch
