import astropy.io.fits
from astropy.io.fits import Header
import numpy as np
from abc import ABC
from collections.abc import Callable

import pandas as pd

from winterdrp.processors.base_processor import BaseImageProcessor, BaseDataframeProcessor
import logging
from winterdrp.processors.database.base_database_processor import BaseDatabaseProcessor, DataBaseError
from winterdrp.processors.database.postgres import import_from_db

logger = logging.getLogger(__name__)


class BaseDatabaseImporter(BaseDatabaseProcessor, ABC):

    base_key = "dbimporter"


def update_header_with_single_match(
        header: Header,
        res: list[dict]
) -> Header:
    assert len(res) == 1

    for key, value in res[0]:
        header[key] = value

    return header


class BaseImageDatabaseImporter(BaseDatabaseImporter, BaseImageProcessor):

    def __init__(
            self,
            db_output_columns: str | list[str],
            output_alias_map: str | list[str] = None,
            update_header: Callable[[Header, list[dict]], Header] = update_header_with_single_match,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.update_header = update_header
        self.db_output_columns = db_output_columns
        self.output_alias_map = output_alias_map

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, header in enumerate(headers):

            query_columns, accepted_values = self.get_constraints(header)

            res = import_from_db(
                db_name=self.db_name,
                db_table=self.db_table,
                db_query_columns=query_columns,
                db_accepted_values=accepted_values,
                db_output_columns=self.db_output_columns,
                output_alias_map=self.output_alias_map,
                db_user=self.db_user,
                password=self.db_password
            )

            headers[i] = self.update_header(header, res)

        return images, headers

    def get_constraints(self, header):
        raise NotImplementedError


class CrossmatchDatabaseWithHeader(BaseImageDatabaseImporter):

    def __init__(
            self,
            db_query_columns: str | list[str],
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.db_query_columns = db_query_columns

    def get_constraints(self, header) -> list[str]:
        accepted_values = [header[x.upper()] for x in self.db_query_columns]
        return accepted_values





class DatabaseDataframeImporter(BaseDatabaseImporter, BaseDataframeProcessor):
    pass

