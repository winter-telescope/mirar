"""
Module to generate a CSV log of observations
"""
import logging
import os
from typing import Optional

import pandas as pd

from mirar.data import ImageBatch
from mirar.paths import BASE_NAME_KEY, core_fields, get_output_path
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

default_log_keys = [BASE_NAME_KEY] + core_fields


class CSVLog(BaseImageProcessor):
    """
    Processor to generate a CSV log
    """

    base_key = "csvlog"

    def __init__(
        self,
        export_keys: Optional[list[str]] = None,
        output_sub_dir: str = "",
        output_base_dir: Optional[str] = None,
    ):
        super().__init__()
        if export_keys is None:
            export_keys = default_log_keys
        self.export_keys = export_keys
        self.output_sub_dir = output_sub_dir
        self.output_base_dir = output_base_dir

    def __str__(self) -> str:
        return "Processor to create a CSV log summarising the image metadata."

    def get_log_name(self) -> str:
        """
        Returns the custom log name

        :return: Lof file name
        """
        return f"{self.night}_log.csv"

    def get_output_path(self) -> str:
        """
        Returns the full log output path

        :return: log path
        """
        output_base_dir = self.output_base_dir
        if output_base_dir is None:
            output_base_dir = self.night_sub_dir

        output_path = get_output_path(
            base_name=self.get_log_name(),
            dir_root=output_base_dir,
            sub_dir=self.output_sub_dir,
        )

        try:
            os.makedirs(os.path.dirname(output_path))
        except OSError:
            pass

        return output_path

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        output_path = self.get_output_path()

        all_rows = []

        for image in batch:
            row = []
            for key in self.export_keys:
                row.append(image[key])

            all_rows.append(row)

        log = pd.DataFrame(all_rows, columns=self.export_keys)

        logger.info(f"Saving log to: {output_path}")
        log.to_csv(output_path)

        return batch
