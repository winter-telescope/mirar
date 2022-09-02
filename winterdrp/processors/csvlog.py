import astropy.io.fits
import numpy as np
import os
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import core_fields, base_name_key, get_output_path
import logging
import pandas as pd

logger = logging.getLogger(__name__)

default_keys = [
    base_name_key
] + core_fields


class CSVLog(BaseImageProcessor):

    base_key = "csvlog"

    def __init__(
            self,
            export_keys: list[str] = default_keys,
            output_sub_dir: str = "",
            output_base_dir: str = None,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.export_keys = export_keys
        self.output_sub_dir = output_sub_dir
        self.output_base_dir = output_base_dir

    def get_log_name(self):
        return f"{self.night}_log.csv"

    def get_output_path(self):
        output_base_dir = self.output_base_dir
        if output_base_dir is None:
            output_base_dir = self.night_sub_dir

        output_path = get_output_path(
            base_name=self.get_log_name(),
            dir_root=output_base_dir,
            sub_dir=self.output_sub_dir
        )

        try:
            os.makedirs(os.path.dirname(output_path))
        except OSError:
            pass

        return output_path

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        output_path = self.get_output_path()

        all_rows = []

        for header in headers:
            row = []
            for key in self.export_keys:
                row.append(header[key])

            all_rows.append(row)

        log = pd.DataFrame(all_rows, columns=self.export_keys)

        logger.info(f"Saving log to: {output_path}")
        log.to_csv(output_path)

        return images, headers
