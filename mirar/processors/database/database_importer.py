"""
Module containing processors which import values from a database
"""
import logging
from abc import ABC
from collections.abc import Callable
from typing import Optional

import pandas as pd

from mirar.data import DataBlock, Image, ImageBatch, SourceBatch
from mirar.processors.base_processor import BaseDataframeProcessor, BaseImageProcessor
from mirar.processors.database.base_database_processor import BaseDatabaseProcessor
from mirar.processors.database.constraints import DBQueryConstraints

logger = logging.getLogger(__name__)


class BaseDatabaseImporter(BaseDatabaseProcessor, ABC):
    """Base Class for any database importer"""

    base_key = "dbimporter"

    def __init__(self, *args, boolean_match_key: Optional[str] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.boolean_match_key = boolean_match_key


def update_header_with_single_match(data: DataBlock, res: list[dict]) -> DataBlock:
    """
    Update a datablock with a single db query result

    :param data: datablock to update
    :param res: corresponding db query
    :return: updated datablock
    """
    assert len(res) == 1

    for key, value in res[0]:
        data[key] = value

    return data


class BaseImageDatabaseImporter(BaseDatabaseImporter, BaseImageProcessor):
    """
    Processor to import data from images
    """

    def __init__(
        self,
        db_output_columns: str | list[str],
        output_alias_map: Optional[str | list[str]] = None,
        update_header: Callable[
            [Image, list[dict]], Image
        ] = update_header_with_single_match,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.update_header = update_header
        self.db_output_columns = db_output_columns
        self.output_alias_map = output_alias_map

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for i, image in enumerate(batch):
            query_constraints = self.get_constraints(image)

            res = self.pg_user.import_from_db(
                db_name=self.db_name,
                db_table=self.db_table,
                db_constraints=query_constraints,
                db_output_columns=self.db_output_columns,
                output_alias_map=self.output_alias_map,
            )

            image = self.update_header(image, res)

            if self.boolean_match_key is not None:
                image[self.boolean_match_key] = len(res) > 0

            batch[i] = image

        return batch

    def get_constraints(self, data: DataBlock) -> None | DBQueryConstraints:
        """
        Get db query constraints for a given datablock object

        :param data: data block
        :return: db query constraints object
        """
        raise NotImplementedError()


class CrossmatchDatabaseWithHeader(BaseImageDatabaseImporter):
    """Processor to crossmatch to a database"""

    def __init__(self, db_query_columns: str | list[str], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.db_query_columns = db_query_columns

    def get_accepted_values(self, data: DataBlock) -> list[str | float | int]:
        """
        Get list of accepted values for crossmatch query

        :param data: datablock
        :return: accepted values from datablock
        """
        accepted_values = [data[x.upper()] for x in self.db_query_columns]
        return accepted_values

    def get_constraints(self, data: DataBlock) -> DBQueryConstraints:
        """
        Get db query constraints for a datablock

        :param data: datablock
        :return: list of constraints
        """
        query_columns = self.db_query_columns
        accepted_values = self.get_accepted_values(data)
        comparison_types = ["=" for _ in self.db_query_columns]
        query_constraints = DBQueryConstraints(
            columns=query_columns,
            accepted_values=accepted_values,
            comparison_types=comparison_types,
        )
        return query_constraints


# def update_dataframe_with_single_match(
#     candidate_table: pd.DataFrame, results: list[dict]
# ) -> pd.DataFrame:
#     for res in results:
#         assert len(res) < 1
#
#     keys = results[0].keys()
#     for key in keys:
#         candidate_table[key] = [x[0][key] for x in results]
#
#     return candidate_table


class DatabaseDataframeImporter(BaseDatabaseImporter, BaseDataframeProcessor, ABC):
    """
    Base Class for dataframe DB importers
    """

    def __init__(
        self,
        db_output_columns: str | list[str],
        output_alias_map: Optional[str | list[str]] = None,
        max_num_results: Optional[int] = None,
        **kwargs,
    ):
        self.db_output_columns = db_output_columns
        self.output_alias_map = output_alias_map
        self.max_num_results = max_num_results
        super().__init__(**kwargs)

    # def _apply_to_candidates(
    #     self,
    #     batch: SourceBatch,
    # ) -> SourceBatch:
    #
    #     for source_table in batch:
    #         candidate_table = source_table.get_data()
    #         results = []
    #         for ind in range(len(candidate_table)):
    #             cand = candidate_table.loc[ind]
    #             (
    #                 query_columns,
    #                 comparison_values,
    #                 comparison_types,
    #             ) = self.get_constraints(cand)
    #             res = import_from_db(
    #                 db_name=self.db_name,
    #                 db_table=self.db_table,
    #                 columns=query_columns,
    #                 accepted_values=comparison_values,
    #                 comparison_types=comparison_types,
    #                 db_output_columns=self.db_output_columns,
    #                 output_alias_map=self.output_alias_map,
    #                 db_user=self.db_user,
    #                 password=self.db_password,
    #                 max_num_results=self.max_num_results,
    #             )
    #             results.append(res)
    #         new_df = self.update_dataframe(candidate_table, results)
    #         source_table.set_data(new_df)
    #     return batch
    #
    # def get_constraints(self, batch: SourceBatch):
    #     raise NotImplementedError


class DatabaseCrossmatchImporter(DatabaseDataframeImporter, BaseDataframeProcessor):
    """
    Processor to crossmatch to sources in a database
    """

    def __init__(
        self,
        crossmatch_radius_arcsec: float,
        user_defined_constraints: Optional[
            Callable[[pd.DataFrame], DBQueryConstraints]
        ] = None,
        ra_field_name: str = "ra",
        dec_field_name: str = "dec",
        order_field_name: Optional[str] = None,
        order_ascending: bool = False,
        q3c_bool: bool = False,
        query_dist: bool = False,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.xmatch_radius_arcsec = crossmatch_radius_arcsec
        self.ra_field_name = ra_field_name
        self.dec_field_name = dec_field_name
        self.q3c = q3c_bool
        self.user_defined_constraints = user_defined_constraints
        self.order_field_name = order_field_name
        self.order_ascending = order_ascending
        self.query_dist = query_dist

    def update_dataframe(
        self, candidate_table: pd.DataFrame, results: list[list[dict]]
    ) -> pd.DataFrame:
        """
        Update a dataframe with db results

        :param candidate_table: pandas table
        :param results: results from db query
        :return: updated dataframe
        """
        assert len(results) == len(candidate_table)
        keys = results[0][0].keys()
        for num in range(len(results[0])):
            for key in keys:
                candidate_table[f"{key}{num + 1}"] = [x[num][key] for x in results]
        return candidate_table

    def get_source_constraints(self, cand: pd.DataFrame) -> None | DBQueryConstraints:
        """
        Get db query constraints for a single source

        :param cand: single source
        :return: db constraint
        """
        if self.user_defined_constraints is None:
            return None
        return self.user_defined_constraints(cand)

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            candidate_table = source_table.get_data()
            results = []
            for ind in range(len(candidate_table)):
                cand = candidate_table.loc[ind]
                query_constraints = self.get_source_constraints(cand)
                res = self.pg_user.crossmatch_with_database(
                    db_name=self.db_name,
                    db_table=self.db_table,
                    db_output_columns=self.db_output_columns,
                    output_alias_map=self.output_alias_map,
                    query_constraints=query_constraints,
                    ra=cand[self.ra_field_name],
                    dec=cand[self.dec_field_name],
                    crossmatch_radius_arcsec=self.xmatch_radius_arcsec,
                    query_distance_bool=self.query_dist,
                    q3c_bool=self.q3c,
                    order_field_name=self.order_field_name,
                    num_limit=self.max_num_results,
                )
                results.append(res)
            new_table = self.update_dataframe(candidate_table, results)
            source_table.set_data(new_table)
        return batch


class DatabaseHistoryImporter(DatabaseCrossmatchImporter):
    """
    Processor to import previous detections of a source from a database
    """

    def __init__(
        self,
        history_duration_days: float,
        time_field_name: str = "jd",
        history_key: str = "prv_candidates",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.history_key = history_key
        self.history_duration_days = history_duration_days
        self.time_field_name = time_field_name
        logger.info(f"Update db is {self.update_dataframe}")

    def update_dataframe(
        self,
        candidate_table: pd.DataFrame,
        results: list[list[dict]],
    ) -> pd.DataFrame:
        """
        Update a pandas dataframe with the number of previous detections

        :param candidate_table: Pandas dataframe
        :param results: db query results
        :return: updated pandas dataframe
        """
        assert len(results) == len(candidate_table)
        candidate_table[self.history_key] = results
        return candidate_table

    def get_source_constraints(self, cand: pd.DataFrame) -> DBQueryConstraints:
        t_detection = float(cand[self.time_field_name])

        query_constraints = DBQueryConstraints(
            columns=self.time_field_name,
            accepted_values=t_detection - self.history_duration_days,
            comparison_types=">",
        )
        if self.user_defined_constraints is not None:
            query_constraints += self.user_defined_constraints(cand)
        return query_constraints
