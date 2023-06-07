"""
Module for processors to modify database entries
"""
import logging
from abc import ABC
from typing import Optional

from mirar.data import ImageBatch
from mirar.processors.database.constraints import DBQueryConstraints
from mirar.processors.database.database_importer import (
    BaseDatabaseImporter,
    BaseImageDatabaseImporter,
)

logger = logging.getLogger(__name__)


class BaseDatabaseModifier(BaseDatabaseImporter, ABC):
    """
    Base class for database modifiers
    """

    base_key = "dbmodifier"

    def __init__(self, db_alter_columns: Optional[str] = None, **kwargs):
        super().__init__(db_output_columns=db_alter_columns, **kwargs)
        self.db_alter_columns = db_alter_columns


class ImageDatabaseModifier(BaseDatabaseModifier, BaseImageDatabaseImporter, ABC):
    """Base Class for modifying image entries in a database"""

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            query_constraints = self.get_constraints(image)

            self.pg_user.modify_db_entry(
                value_dict=image,
                db_constraints=query_constraints,
                db_alter_columns=self.db_alter_columns,
                db_table=self.db_table,
                db_name=self.db_name,
            )

        return batch


class ModifyImageDatabaseSeq(ImageDatabaseModifier):
    """Processor to modify images in a database with a sequence"""

    def __init__(self, sequence_key: Optional[str | list[str]] = None, **kwargs):
        super().__init__(**kwargs)
        self.sequence_key = sequence_key

    def get_constraints(self, data):
        if self.sequence_key is None:
            self.sequence_key = list(
                self.pg_user.get_sequence_keys_from_table(self.db_table, self.db_name)
            )

        accepted_values = [data[x.upper()] for x in self.sequence_key]
        comparison_types = ["="] * len(accepted_values)

        query_constraints = DBQueryConstraints(
            columns=self.sequence_key,
            accepted_values=accepted_values,
            comparison_types=comparison_types,
        )

        return query_constraints
