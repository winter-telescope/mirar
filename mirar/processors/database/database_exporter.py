"""
Module containing processors for exporting to databases
"""
import logging
from abc import ABC

from mirar.data import ImageBatch, SourceBatch
from mirar.processors.base_processor import BaseDataframeProcessor, BaseImageProcessor
from mirar.processors.database.base_database_processor import BaseDatabaseProcessor

logger = logging.getLogger(__name__)


class BaseDatabaseExporter(BaseDatabaseProcessor, ABC):
    """
    Base class for DB exporters
    """

    base_key = "dbexporter"

    def __str__(self):
        return (
            f"Processor to save "
            f"{['candidates', 'images'][isinstance(self, BaseImageProcessor)]} "
            f"to the '{self.db_table}' table of the '{self.db_name}' Postgres database."
        )


class DatabaseImageExporter(BaseDatabaseExporter, BaseImageProcessor):
    """
    Processor for exporting images to a database
    """

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image in batch:
            primary_keys, primary_key_values = self.pg_user.export_to_db(
                image,
                db_name=self.db_name,
                db_table=self.db_table,
                duplicate_protocol=self.duplicate_protocol,
            )

            for ind, key in enumerate(primary_keys):
                image[key] = primary_key_values[ind]
        return batch


class DatabaseDataframeExporter(BaseDatabaseExporter, BaseDataframeProcessor):
    """
    Processor for exporting sources to a database
    """

    def _apply_to_candidates(self, batch: SourceBatch) -> SourceBatch:
        for source_list in batch:
            candidate_table = source_list.get_data()

            primary_key_dict = {}
            for _, candidate_row in candidate_table.iterrows():
                primary_keys, primary_key_values = self.pg_user.export_to_db(
                    candidate_row.to_dict(),
                    db_name=self.db_name,
                    db_table=self.db_table,
                    duplicate_protocol=self.duplicate_protocol,
                )
                for ind, key in enumerate(primary_keys):
                    if key not in primary_key_dict:
                        primary_key_dict[key] = [primary_key_values[ind]]
                    else:
                        primary_key_dict[key].append(primary_key_values[ind])

            for key, val in primary_key_dict.items():
                candidate_table[key] = val

            source_list.set_data(candidate_table)

        return batch
