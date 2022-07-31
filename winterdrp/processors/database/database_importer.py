import astropy.io.fits
from astropy.io.fits import Header
import numpy as np
from abc import ABC
from collections.abc import Callable

import pandas as pd

from winterdrp.processors.base_processor import BaseImageProcessor, BaseDataframeProcessor
import logging
from winterdrp.processors.database.base_database_processor import BaseDatabaseProcessor, DataBaseError
from winterdrp.processors.database.postgres import import_from_db, xmatch_import_db

logger = logging.getLogger(__name__)


class BaseDatabaseImporter(BaseDatabaseProcessor, ABC):
    base_key = "dbimporter"

    def __init__(
            self,
            boolean_match_key: str = None,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.boolean_match_key = boolean_match_key


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

            new_header = self.update_header(header, res)

            if self.boolean_match_key is not None:
                new_header[self.boolean_match_key] = len(res) > 0

            headers[i] = new_header

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


def update_dataframe_with_single_match(
        candidate_table: pd.DataFrame,
        results: list[dict]
) -> pd.DataFrame:
    for res in results:
        assert len(res) < 1

    keys = results[0].keys()
    for key in keys:
        candidate_table[key] = [x[0][key] for x in results]

    return candidate_table


class DatabaseDataframeImporter(BaseDatabaseImporter, BaseDataframeProcessor):
    def __init__(self,
                 db_output_columns: str | list[str],
                 output_alias_map: str | list[str] = None,
                 update_dataframe: Callable[
                     [pd.DataFrame, list[list[dict]]], pd.DataFrame] = update_dataframe_with_single_match,
                 max_num_results: int = None,
                 *args, **kwargs):
        self.db_output_columns = db_output_columns
        self.output_alias_map = output_alias_map
        self.update_dataframe = update_dataframe
        self.max_num_results = max_num_results
        super(DatabaseDataframeImporter, self).__init__(*args, **kwargs)

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        results = []
        for ind in range(len(candidate_table)):
            cand = candidate_table.loc[ind]
            query_columns, comparison_values, comparison_types = self.get_constraints(cand)
            res = import_from_db(
                db_name=self.db_name,
                db_table=self.db_table,
                db_query_columns=query_columns,
                db_accepted_values=comparison_values,
                db_output_columns=self.db_output_columns,
                output_alias_map=self.output_alias_map,
                db_user=self.db_user,
                password=self.db_password,
                max_num_results=self.max_num_results
            )

            results.append(res)
        new_df = self.update_dataframe(candidate_table, results)
        return new_df

    def get_constraints(self, cand):
        raise NotImplementedError


def no_additional_constraints(cand):
    return [], [], []


def update_xmatch_dataframe(
        candidate_table: pd.DataFrame,
        results: list[list[dict]]) -> pd.DataFrame:
    assert len(results) == len(candidate_table)
    keys = results[0][0].keys()
    for num in range(len(results[0])):
        for key in keys:
            candidate_table[f"{key}{num + 1}"] = [x[num][key] for x in results]
    return candidate_table


class DatabaseXMatchImporter(DatabaseDataframeImporter, BaseDataframeProcessor):
    def __init__(self,
                 xmatch_radius_arcsec: float,
                 user_defined_constraints: Callable[[pd.DataFrame], tuple] = no_additional_constraints,
                 update_dataframe: Callable[
                     [pd.DataFrame, list[list[dict]]], pd.DataFrame] = update_xmatch_dataframe,
                 ra_field_name: str = 'ra',
                 dec_field_name: str = 'dec',
                 order_field_name: str = None,
                 order_ascending: bool = False,
                 q3c: bool = False,
                 query_dist: bool = False,
                 *args, **kwargs):
        super(DatabaseXMatchImporter, self).__init__(*args, **kwargs)
        self.xmatch_radius_arcsec = xmatch_radius_arcsec
        self.ra_field_name = ra_field_name
        self.dec_field_name = dec_field_name
        self.q3c = q3c
        self.user_defined_constraints = user_defined_constraints
        self.order_field_name = order_field_name
        self.order_ascending = order_ascending
        self.query_dist = query_dist
        self.update_dataframe = update_dataframe

    def get_constraints(self, cand):
        query_columns, comparison_types, accepted_values = self.user_defined_constraints(cand)
        return query_columns, comparison_types, accepted_values

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        results = []
        for ind in range(len(candidate_table)):
            cand = candidate_table.loc[ind]
            query_columns, comparison_types, accepted_values = self.get_constraints(cand)
            res = xmatch_import_db(db_name=self.db_name,
                                   db_table=self.db_table,
                                   db_user=self.db_user,
                                   db_password=self.db_password,
                                   db_output_columns=self.db_output_columns,
                                   output_alias_map=self.output_alias_map,
                                   db_query_columns=query_columns,
                                   db_comparison_types=comparison_types,
                                   db_accepted_values=accepted_values,
                                   ra=cand[self.ra_field_name],
                                   dec=cand[self.dec_field_name],
                                   xmatch_radius_arcsec=self.xmatch_radius_arcsec,
                                   query_dist=self.query_dist,
                                   q3c=self.q3c,
                                   order_field_name=self.order_field_name,
                                   order_ascending=self.order_ascending,
                                   num_limit=self.max_num_results,
                                   )
            results.append(res)

        new_table = self.update_dataframe(candidate_table, results)
        return new_table


def update_history_dataframe(
        candidate_table: pd.DataFrame,
        results: list[list[dict]]
):
    assert len(results) == len(candidate_table)
    candidate_table['prv_candidates'] = results
    return candidate_table


class DatabaseHistoryImporter(DatabaseXMatchImporter):
    def __init__(self,
                 history_duration_days: float,
                 time_field_name: str = 'jd',
                 update_dataframe: Callable[
                     [pd.DataFrame, list[list[dict]]], pd.DataFrame] = update_history_dataframe,
                 *args, **kwargs):
        super(DatabaseHistoryImporter, self).__init__(*args, **kwargs)
        self.history_duration_days = history_duration_days
        self.time_field_name = time_field_name
        self.update_dataframe = update_dataframe
        logger.info(f'Update db is {self.update_dataframe}')

    def get_constraints(self, cand):
        query_columns, comparison_types, accepted_values = self.user_defined_constraints(cand)
        query_columns.append(self.time_field_name)
        comparison_types.append('>')
        accepted_values.append(cand[self.time_field_name] - self.history_duration_days)
        return query_columns, comparison_types, accepted_values
