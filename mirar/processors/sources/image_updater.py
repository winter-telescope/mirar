"""
Module to update the header of images with data in a source table
"""

import json
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from mirar.data import SourceBatch
from mirar.io import open_raw_image, save_fits
from mirar.paths import BASE_NAME_KEY, base_output_dir, get_output_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class ImageUpdater(BaseSourceProcessor):
    """
    Update the header of images with data in a source table
    """

    base_key = "IMGUPDATE"

    def __init__(
        self,
        modify_dir_name: Optional[str] = None,
        input_dir: str | Path = base_output_dir,
        update_keys: list[str] | str | None = None,
        overwrite: bool = False,
        include_table: bool = True,
    ):
        super().__init__()
        self.modify_dir_name = modify_dir_name
        self.input_dir = Path(input_dir)

        if update_keys is not None:
            if not isinstance(update_keys, list):
                update_keys = [update_keys]

        self.update_keys = update_keys
        self.overwrite = overwrite
        self.include_table = include_table

    def description(self) -> str:
        return (
            f"Processor to update image headers in '{self.modify_dir_name}' directory."
        )

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        mod_dir = get_output_dir(
            dir_root=self.modify_dir_name,
            sub_dir=self.night_sub_dir,
            output_dir=self.input_dir,
        )

        for source_list in batch:

            image_path = mod_dir / source_list[BASE_NAME_KEY]

            source_table = source_list.get_data()
            metadata = source_list.get_metadata()

            if self.include_table:
                if len(source_table) != 1:
                    err = (
                        f"Can only include source table if it has a single row, "
                        f"but this table has {len(source_table)} rows."
                    )
                    logger.error(err)
                    raise ValueError(err)

                source_row = source_table.iloc[0]
                export_dict = self.generate_super_dict(metadata, source_row)
            else:
                export_dict = metadata

            if self.update_keys is not None:
                export_dict = {k: export_dict[k] for k in self.update_keys}

            export_dict = json.loads(pd.Series(export_dict).to_json())

            img = open_raw_image(image_path)
            header = img.get_header()

            for key, value in export_dict.items():
                if (key in header) & (not self.overwrite):
                    logger.debug(f"Skipping overwrite of {key} in {image_path}")
                    continue

                if (key not in header) | self.overwrite:
                    try:
                        header[key] = value
                    except ValueError:
                        logger.debug(f"Failed to copy {key} in {image_path}")
                        continue

            img.set_header(header)

            logger.debug(f"Updating header of {image_path}")

            save_fits(img, image_path)

        return batch
