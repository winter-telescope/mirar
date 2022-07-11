import astropy.io.fits
import numpy as np
from abc import ABC

import pandas as pd

from winterdrp.processors.base_processor import BaseImageProcessor, BaseDataframeProcessor
import logging
from winterdrp.processors.database.postgres import export_to_db
from winterdrp.processors.database.base_database_processor import BaseDatabaseProcessor

logger = logging.getLogger(__name__)


class BaseDatabaseImporter(BaseDatabaseProcessor, ABC):

    base_key = "dbimporter"

    def __init__(
            self,
            result_field_name: str,
            import_from_db: callable,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.result_field_name = result_field_name
        self.import_from_db = import_from_db


class DatabaseImageImporter(BaseDatabaseImporter, BaseImageProcessor):

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, header in enumerate(headers):
            header = self.import_from_db(
                header,
                db_name=self.db_name,
                db_table=self.db_table,
                db_user=self.db_user,
                password=self.db_password
            )
            headers[i] = header

        return images, headers


class DatabaseDataframeImporter(BaseDatabaseImporter, BaseDataframeProcessor):

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame
    ) -> pd.DataFrame:

        new_table = pd.DataFrame()

        for index, candidate_row in candidate_table.iterrows():
            new_row = self.import_from_db(
                candidate_row,
                db_name=self.db_name,
                db_table=self.db_table,
                db_user=self.db_user,
                password=self.db_password
            )

            new_table.append(new_row, ignore_index=True)

        return new_table

