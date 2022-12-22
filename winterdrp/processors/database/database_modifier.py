"""
Module for processors to modify database entries
"""
import logging
from abc import ABC
from typing import Optional

from winterdrp.data import ImageBatch
from winterdrp.processors.database.database_importer import (
    BaseDatabaseImporter,
    BaseImageDatabaseImporter,
)
from winterdrp.processors.database.postgres import (
    get_sequence_keys_from_table,
    modify_db_entry,
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
            query_columns, accepted_values, accepted_types = self.get_constraints(image)
            logger.debug(f"{query_columns}, {accepted_values}, {accepted_types}")

            modify_db_entry(
                value_dict=image,
                db_query_columns=query_columns,
                db_query_values=accepted_values,
                db_query_comparison_types=accepted_types,
                db_alter_columns=self.db_alter_columns,
                db_table=self.db_table,
                db_name=self.db_name,
                db_user=self.db_user,
                password=self.db_password,
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
                get_sequence_keys_from_table(
                    self.db_table, self.db_name, self.db_user, self.db_password
                )
            )

        accepted_values = [data[x.upper()] for x in self.sequence_key]
        accepted_types = ["="] * len(accepted_values)
        return self.sequence_key, accepted_values, accepted_types
