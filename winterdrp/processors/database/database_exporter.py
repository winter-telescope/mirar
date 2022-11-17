import astropy.io.fits
import numpy as np
from abc import ABC

import pandas as pd

from winterdrp.processors.base_processor import BaseImageProcessor, BaseDataframeProcessor
import logging
from winterdrp.processors.database.postgres import DataBaseError, export_to_db
from winterdrp.processors.database.base_database_processor import BaseDatabaseProcessor
from winterdrp.data import ImageBatch


logger = logging.getLogger(__name__)


class BaseDatabaseExporter(BaseDatabaseProcessor, ABC):
    base_key = "dbexporter"

    def __str__(self):
        return f"Processor to save {['candidates', 'images'][isinstance(self, BaseImageProcessor)]} " \
               f"to the '{self.db_table}' table of the '{self.db_name}' Postgres database."


class DatabaseImageExporter(BaseDatabaseExporter, BaseImageProcessor):

    def _apply_to_images(
            self,
            batch: ImageBatch
    ) -> ImageBatch:

        for image in batch:
            primary_keys, primary_key_values = export_to_db(
                image,
                db_name=self.db_name,
                db_table=self.db_table,
                db_user=self.db_user,
                password=self.db_password,
                duplicate_protocol=self.duplicate_protocol
            )

            for ind, key in enumerate(primary_keys):
                image[key] = primary_key_values[ind]
        return batch


class DatabaseDataframeExporter(BaseDatabaseExporter, BaseDataframeProcessor):

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame
    ) -> pd.DataFrame:

        primary_key_dict = {}
        for index, candidate_row in candidate_table.iterrows():
            primary_keys, primary_key_values = export_to_db(
                candidate_row.to_dict(),
                db_name=self.db_name,
                db_table=self.db_table,
                db_user=self.db_user,
                password=self.db_password
            )
            for ind, key in enumerate(primary_keys):
                if key not in primary_key_dict:
                    primary_key_dict[key] = [primary_key_values[ind]]
                else:
                    primary_key_dict[key].append(primary_key_values[ind])
                # candidate_row[key] = primary_key_values[ind]

            # new_table = pd.concat([new_table,candidate_row])
            # new_table.append(candidate_row, ignore_index=True)
        for k in primary_key_dict:
            candidate_table[k] = primary_key_dict[k]

        return candidate_table
