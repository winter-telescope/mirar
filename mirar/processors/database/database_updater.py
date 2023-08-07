"""
Module for processors to modify database entries
"""
import logging
from abc import ABC
from typing import Optional

import pandas as pd

from mirar.data import ImageBatch
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions import select_from_table
from mirar.database.utils import get_sequence_key_names_from_table
from mirar.processors.database.database_inserter import DatabaseImageInserter
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
            val_dict = self.generate_value_dict(image)
            new = self.db_table(**val_dict)
            new.update_entry(update_keys=self.db_alter_columns)
        return batch

    @staticmethod
    def generate_value_dict(image):
        """
        Get the value dictionary for an image

        :param image: Image
        :return: Value dictionary
        """
        return DatabaseImageInserter.generate_value_dict(image)


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

        :param data: Image
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


class ImageDatabaseMultiEntryUpdater(ImageSequenceDatabaseUpdater):
    """
    Processor to modify multiple entries specified by a list of sequences in an
    image database
    """

    def __init__(self, sequence_key: str, **kwargs):
        super().__init__(**kwargs)
        self.sequence_key = sequence_key.lower()

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            try:
                unique_key_vals = [int(y) for y in image[self.sequence_key].split(",")]
            except ValueError as exc:
                raise ValueError("Sequence keys must be integers") from exc

            data_df = pd.DataFrame(unique_key_vals, columns=[self.sequence_key])

            for _, row in data_df.iterrows():
                constraints = DBQueryConstraints(
                    columns=self.sequence_key, accepted_values=row[self.sequence_key]
                )

                old = select_from_table(
                    db_constraints=constraints,
                    sql_table=self.db_table.sql_model,
                )

                assert (
                    len(old) == 1
                ), f"Multiple entries found for unique key {self.sequence_key}"

                for key in self.db_alter_columns:
                    old[key] = image[key]

                new = self.db_table(**old.to_dict(orient="records")[0])

                new.update_entry(update_keys=self.db_alter_columns)

        return batch
