"""
Module for saving images
"""
import shutil
from pathlib import Path

from mirar.data import ImageBatch
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_SAVE_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    base_output_dir,
    get_output_path,
)
from mirar.processors.base_processor import BaseImageProcessor


class ImageSaver(BaseImageProcessor):
    """
    Processor to save images
    """

    base_key = "save"

    def __init__(
        self,
        output_dir_name: str,
        write_mask: bool = True,
        output_dir: str | Path = base_output_dir,
    ):
        super().__init__()
        self.output_dir_name = output_dir_name
        self.write_mask = write_mask
        self.output_dir = Path(output_dir)

    def __str__(self):
        return f"Processor to save images to the '{self.output_dir_name}' subdirectory"

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            path = get_output_path(
                image[BASE_NAME_KEY],
                dir_root=self.output_dir_name,
                sub_dir=self.night_sub_dir,
                output_dir=self.output_dir,
            )

            path.parent.mkdir(parents=True, exist_ok=True)

            image[LATEST_SAVE_KEY] = str(path)
            if self.write_mask:
                weight_image_found, mask_path = False, ""
                if LATEST_WEIGHT_SAVE_KEY in image.header.keys():
                    existing_weightpath = Path(image[LATEST_WEIGHT_SAVE_KEY])
                    if existing_weightpath.exists():
                        weight_image_found = True
                        mask_path = get_output_path(
                            existing_weightpath.name,
                            dir_root=self.output_dir_name,
                            sub_dir=self.night_sub_dir,
                            output_dir=self.output_dir,
                        )
                        if existing_weightpath != mask_path:
                            shutil.copy(existing_weightpath, mask_path)

                if not weight_image_found:
                    mask_path = self.save_weight_image(image, img_path=path)
                image[LATEST_WEIGHT_SAVE_KEY] = str(mask_path)
            self.save_fits(image, path)

        return batch
