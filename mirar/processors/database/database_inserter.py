"""
Module containing processors for exporting to databases
"""

import logging
from abc import ABC

import numpy as np
import pandas as pd

from mirar.data import Image, ImageBatch, SourceBatch
from mirar.errors.exceptions import BaseProcessorError
from mirar.processors.base_processor import (
    BaseImageProcessor,
    BaseSourceProcessor,
    PrerequisiteError,
)
from mirar.processors.database.base_database_processor import BaseDatabaseProcessor
from mirar.processors.utils.image_selector import ImageBatcher

logger = logging.getLogger(__name__)


class ImageBatchDatabaseExporterError(BaseProcessorError):
    """
    Error for ImageBatchDatabaseExporter
    """


class BaseDatabaseInserter(BaseDatabaseProcessor, ABC):
    """
    Base class for DB inserter
    """

    base_key = "dbinserter"
    max_n_cpu = 1

    def __str__(self):
        return (
            f"Processor to save "
            f"{['candidates', 'images'][isinstance(self, BaseImageProcessor)]} "
            f"to the '{self.db_table.__name__}' table of "
            f"the '{self.db_name}' Postgres database."
        )


class DatabaseImageInserter(BaseDatabaseInserter, BaseImageProcessor):
    """
    Processor for exporting images to a database
    """

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image in batch:
            val_dict = self.generate_value_dict(image)

            new = self.db_table(**val_dict)
            res = new.insert_entry()

            assert len(res) == 1

            for key in res.columns:
                image[key] = res[key].iloc[0]
        return batch

    @staticmethod
    def generate_value_dict(image: Image) -> dict:
        """
        Generate a dictionary of lower case keys and values from an image

        :param image: Image
        :return: Dictionary of lower case keys and values
        """
        return {key.lower(): image[key] for key in image.keys()}


class DatabaseSourceInserter(BaseDatabaseInserter, BaseSourceProcessor):
    """
    Processor for exporting sources to a database
    """

    def _apply_to_sources(self, batch: SourceBatch) -> SourceBatch:
        for source_list in batch:
            source_table = source_list.get_data()
            metadata = source_list.get_metadata()

            primary_key_df_list = []
            for _, source_row in source_table.iterrows():
                super_dict = self.generate_super_dict(metadata, source_row)

                new = self.db_table(**super_dict)
                res = new.insert_entry()

                assert len(res) == 1

                primary_key_df_list.append(res.loc[0])

            primary_key_df = pd.DataFrame(primary_key_df_list).reset_index(drop=True)
            for key in primary_key_df:
                source_table[key] = primary_key_df[key]

            source_list.set_data(source_table)

        return batch


class DatabaseImageBatchInserter(DatabaseImageInserter):
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

        val_dict = {key.lower(): image[key] for key in image.keys()}

        new = self.db_table(**val_dict)
        res = new.insert_entry()

        assert len(res) == 1

        for key in res.columns:
            val = res[key].iloc[0]
            for image in batch:
                image[key] = val

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
            raise PrerequisiteError(err)
