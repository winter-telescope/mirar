"""
Module containing processors for exporting to databases
"""
import logging
from abc import ABC

import numpy as np

from mirar.data import ImageBatch, SourceBatch
from mirar.errors.exceptions import BaseProcessorError
from mirar.processors.base_processor import BaseDataframeProcessor, BaseImageProcessor
from mirar.processors.sqldatabase.base_database_processor import BaseDatabaseProcessor
from mirar.processors.utils.image_selector import ImageBatcher

logger = logging.getLogger(__name__)


class ImageBatchDatabaseExporterError(BaseProcessorError):
    """
    Error for ImageBatchDatabaseExporter
    """


class BaseDatabaseExporter(BaseDatabaseProcessor, ABC):
    """
    Base class for DB exporters
    """

    base_key = "dbexporter"
    max_n_cpu = 1

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
                    db_table=self.db_table,
                )
                for ind, key in enumerate(primary_keys):
                    if key not in primary_key_dict:
                        primary_key_dict[key] = [primary_key_values[ind]]
                    else:
                        primary_key_dict[key].append(primary_key_values[ind])
                    # candidate_row[key] = primary_key_values[ind]

                # new_table = pd.concat([new_table,candidate_row])
                # new_table.append(candidate_row, ignore_index=True)
            for key, val in primary_key_dict.items():
                candidate_table[key] = val

            source_list.set_data(candidate_table)

        return batch


class DatabaseImageBatchExporter(DatabaseImageExporter):
    """
    Processor for creating a single entry per batch of images in a database
    """

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        column_names = [
            x for x in self.db_table.__dict__["__annotations__"] if x != "sql_model"
        ]

        for column in column_names:
            try:
                values = [x[column] for x in batch]
            except KeyError as exc:
                err = (
                    f"Key {column} not found in the batch, cannot export it "
                    f"to database. Available keys are {list(batch[0].keys())}"
                )
                logger.error(err)
                raise ImageBatchDatabaseExporterError(err) from exc

            if len(np.unique(values)) > 1:
                err = (
                    f"Key {column} differs across images in the batch, cannot export"
                    f"it to database."
                )
                logger.error(err)
                raise ImageBatchDatabaseExporterError(err)

        image = batch[0]
        logger.debug(f"Trying to export {[image[x] for x in column_names]}")
        primary_keys, primary_key_values = self.pg_user.export_to_db(
            image,
            db_table=self.db_table,
            duplicate_protocol=self.duplicate_protocol,
        )

        for ind, key in enumerate(primary_keys):
            for image in batch:
                image[key] = primary_key_values[ind]
        return batch

    def check_prerequisites(
        self,
    ):
        check = np.sum([isinstance(x, ImageBatcher) for x in self.preceding_steps[-1:]])
        if check < 1:
            err = (
                f"{self.__module__} requires {ImageBatcher} to be run right before it"
                f"as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise ValueError
