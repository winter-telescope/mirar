import astropy.io.fits
from astropy.io.fits import Header
import numpy as np
from abc import ABC
from collections.abc import Callable

import pandas as pd

from winterdrp.processors.base_processor import BaseImageProcessor, BaseDataframeProcessor
import logging
from winterdrp.processors.database.base_database_processor import BaseDatabaseProcessor, DataBaseError
from winterdrp.processors.database.postgres import modify_db_entry, get_sequence_keys_from_table
from winterdrp.processors.database.database_importer import BaseDatabaseImporter, BaseImageDatabaseImporter

logger = logging.getLogger(__name__)


class BaseDatabaseModifier(BaseDatabaseImporter):
    base_key = "dbmodifier"

    def __init__(self,
                 db_alter_columns: str = None,
                 *args,
                 **kwargs):
        super(BaseDatabaseModifier, self).__init__(db_output_columns=db_alter_columns, *args, **kwargs)
        self.db_alter_columns = db_alter_columns


class ImageDatabaseModifier(BaseDatabaseModifier, BaseImageDatabaseImporter):

    def __init__(self,
                 *args,
                 **kwargs):
        super(ImageDatabaseModifier, self).__init__(*args, **kwargs)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:
        for header in headers:
            query_columns, accepted_values, accepted_types = self.get_constraints(header)
            logger.info(f"{query_columns}, {accepted_values}, {accepted_types}")

            query_results = modify_db_entry(value_dict=header,
                                            db_query_columns=query_columns,
                                            db_query_values=accepted_values,
                                            db_query_comparison_types=accepted_types,
                                            db_alter_columns=self.db_alter_columns,
                                            db_table=self.db_table,
                                            db_name=self.db_name,
                                            db_user=self.db_user,
                                            password=self.db_password
                                            )

        return images, headers


class ModifyImageDatabaseSeq(ImageDatabaseModifier):
    def __init__(self,
                 sequence_key: str | list[str] = None,
                 *args,
                 **kwargs):
        super(ModifyImageDatabaseSeq, self).__init__(*args, **kwargs)
        self.sequence_key = sequence_key

    def get_constraints(self, header):
        if self.sequence_key is None:
            self.sequence_key = [x for x in get_sequence_keys_from_table(self.db_table,
                                                                         self.db_name,
                                                                         self.db_user,
                                                                         self.db_password)]

        accepted_values = [header[x.upper()] for x in self.sequence_key]
        accepted_types = ['='] * len(accepted_values)
        return self.sequence_key, accepted_values, accepted_types


class DataframeDatabaseModifier(BaseDatabaseModifier, BaseDataframeProcessor):

    def __init__(self,
                 *args,
                 **kwargs):
        super(DataframeDatabaseModifier, self).__init__(*args, **kwargs)

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        return candidate_table
