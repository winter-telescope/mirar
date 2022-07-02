import astropy.io.fits
import numpy as np
import os
from winterdrp.processors.base_processor import BaseImageProcessor
import logging
from winterdrp.processors.database.postgres import DataBaseError, export_to_db
from winterdrp.processors.database.base_database_processor import BaseDatabaseProcessor

logger = logging.getLogger(__name__)

class DatabaseExporter(BaseDatabaseProcessor, BaseImageProcessor):
    base_key = "dbexport"

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for header in headers:
            primary_keys, primary_key_values = export_to_db(
                header,
                db_name=self.db_name,
                db_table=self.db_table,
                db_user=self.db_user,
                password=self.db_password
            )

            for ind, key in enumerate(primary_keys):
                header[key] = primary_key_values[ind]
        return images, headers

