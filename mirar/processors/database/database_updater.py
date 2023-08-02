"""
Module for processors to modify database entries
"""
import logging
from abc import ABC
from typing import Optional

import numpy as np
import pandas as pd

from mirar.data import Image, ImageBatch
from mirar.database.constraints import DBQueryConstraints
from mirar.database.postgres_utils import get_sequence_key_names_from_table
from mirar.database.transactions import update_database_entry
from mirar.processors.database.database_selector import (
    BaseDatabaseSelector,
    BaseImageDatabaseSelector,
)

logger = logging.getLogger(__name__)


class BaseDatabaseUpdater(BaseDatabaseSelector, ABC):
    """
    Base class for database updaters
    """

    base_key = "dbupdater"

    def __init__(self, db_alter_columns: str | list[str], **kwargs):
        super().__init__(db_output_columns=db_alter_columns, **kwargs)
        if not isinstance(db_alter_columns, list):
            db_alter_columns = [db_alter_columns]
        self.db_alter_columns = db_alter_columns


class ImageDatabaseUpdater(BaseDatabaseUpdater, BaseImageDatabaseSelector, ABC):
    """Base Class for updating image entries in a database"""

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            query_constraints = self.get_constraints(image)

            update_dict = self.get_export_dict(image)

            update_database_entry(
                update_dict=update_dict,
                db_constraints=query_constraints,
                db_model=self.db_table,
            )

        return batch

    def get_export_dict(self, image: Image) -> dict:
        """
        Create a dictionary to export the relevant fields to a database from an image

        :param image: Image to export
        :return: Dictionary of keys/values to export
        """
        new = {}
        for key in self.db_alter_columns:
            new[key] = image[key]
        return new


class ImageSequenceDatabaseUpdater(ImageDatabaseUpdater):
    """
    Processor to modify images in a database with a sequence
    """

    def __init__(self, sequence_key: Optional[str | list[str]] = None, **kwargs):
        super().__init__(**kwargs)
        self.sequence_key = sequence_key

    def get_constraints(self, data) -> DBQueryConstraints:
        """
        Function to get the constraints for a database query

        :param data: Image to get constraints for
        :return: Constraints for a database query
        """
        if self.sequence_key is None:
            self.sequence_key = get_sequence_key_names_from_table(
                self.db_table.sql_model.__tablename__, self.db_name
            )

        accepted_values = [data[x.lower()] for x in self.sequence_key]
        comparison_types = ["="] * len(accepted_values)

        query_constraints = DBQueryConstraints(
            columns=self.sequence_key,
            accepted_values=accepted_values,
            comparison_types=comparison_types,
        )

        return query_constraints


class ImageSequenceDatabaseUpdaterList(ImageSequenceDatabaseUpdater):
    """
    Processor to modify multiple entries specified by a list of sequences in an
    image database
    """

    def __init__(self, sequence_key: Optional[str | list[str]], **kwargs):
        super().__init__(**kwargs)
        self.sequence_key = sequence_key
        if isinstance(self.sequence_key, str):
            self.sequence_key = [self.sequence_key]
        if self.sequence_key is None:
            self.sequence_key = get_sequence_key_names_from_table(
                self.db_table.sql_model.__tablename__, self.db_name
            )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            try:
                data_array = [
                    [int(y) for y in image[x.lower()].split(",")]
                    for x in self.sequence_key
                ]
            except ValueError as exc:
                raise ValueError("Sequence keys must be integers") from exc

            data_array_lens = [len(x) for x in data_array]
            if len(set(data_array_lens)) > 1:
                raise ValueError("All sequence keys must have the same length")

            data_df = pd.DataFrame(
                np.array(data_array).T, columns=[x.lower() for x in self.sequence_key]
            )
            logger.debug(data_df)

            update_dict = self.get_export_dict(image)

            for ind in range(len(data_df)):
                row = data_df.iloc[ind]
                query_constraints = self.get_constraints(row)

                update_database_entry(
                    update_dict=update_dict,
                    db_constraints=query_constraints,
                    db_model=self.db_table,
                )

        return batch
