"""
Module containing postgres util functions
"""
import logging
from typing import Type

import psycopg
from psycopg import errors
from psycopg.rows import Row
from pydantic import ValidationError
from sqlalchemy import inspect, text
from sqlalchemy.exc import IntegrityError
from sqlalchemy_utils import create_database, database_exists

from mirar.data import DataBlock
from mirar.database.base_model import BaseDB
from mirar.database.constants import POSTGRES_DUPLICATE_PROTOCOLS
from mirar.database.engine import get_engine
from mirar.database.user.credentials import (
    DB_PASSWORD,
    DB_PASSWORD_KEY,
    DB_USER,
    DB_USER_KEY,
)
from mirar.errors import ProcessorError

logger = logging.getLogger(__name__)


class DataBaseError(ProcessorError):
    """Error relating to postgres interactions"""


class PostgresUser:
    """
    Basic Postgres user class for executing functions
    """

    user_env_varaiable = DB_USER_KEY
    pass_env_variable = DB_PASSWORD_KEY

    def __init__(self, db_user: str = DB_USER, db_password: str = DB_PASSWORD):
        self.db_user = db_user
        self.db_password = db_password

    def validate_credentials(self):
        """
        Checks that user credentials exist

        :return: None
        """
        if self.db_user is None:
            err = (
                f"'db_user' is set as None. Please pass a db_user as an argument, "
                f"or set the environment variable '{self.user_env_varaiable}'."
            )
            logger.error(err)
            raise DataBaseError(err)

        if self.db_password is None:
            err = (
                f"'db_password' is set as None. Please pass a password as an argument, "
                f"or set the environment variable '{self.pass_env_variable}'."
            )
            logger.error(err)
            raise DataBaseError(err)

        engine = get_engine(
            db_name="postgres",
            db_user=self.db_user,
            db_password=self.db_password,
        )

        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))

    def export_to_db(
        self,
        value_dict: dict | DataBlock,
        db_table: Type[BaseDB],
        duplicate_protocol: str = "fail",
    ) -> tuple[list, list]:
        """
        Export a list of fields in value dict to a batabase table

        :param value_dict: dictionary/DataBlock/other dictonary-like object to export
        :param db_table: table of DB to export to
        :param duplicate_protocol: protocol for handling duplicates,
            in "fail"/"ignore"/"replace"
        :return:
        """

        assert duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

        column_names = [
            x for x in db_table.__dict__["__annotations__"] if x != "sql_model"
        ]

        column_dict = {}
        for column in column_names:
            column_dict[column] = value_dict[column]

        try:
            new = db_table(**column_dict)
        except ValidationError as err:
            logger.error(err)
            raise DataBaseError from err

        db_name = new.sql_model.db_name
        primary_key = inspect(db_table.sql_model).primary_key[0]

        sequence_key_names, sequence_values = [], []
        try:
            sequence_key_names, sequence_values = new.insert_entry()

        except IntegrityError as exc:
            if not isinstance(exc.orig, errors.UniqueViolation):
                raise exc

            if duplicate_protocol == "fail":
                err = (
                    f"Duplicate error, entry with {column_dict} "
                    f"already exists in {db_name}."
                )
                logger.error(err)
                raise errors.UniqueViolation from exc

            if duplicate_protocol == "ignore":
                logger.debug(
                    f"Found duplicate entry in {db_name} - "
                    f"{str(exc)}."
                    f"Ignoring, no new entry made."
                )
                primary_key_val = value_dict[primary_key.name]
                sequence_keys = new.get_sequence_keys()
                sequence_key_names = [k.name for k in sequence_keys]
                sequence_values = []
                if len(sequence_keys) > 0:
                    ret = new.sql_model().select_query(
                        compare_values=[primary_key_val],
                        compare_keys=[primary_key.name],
                        select_keys=sequence_key_names,
                    )
                    sequence_values = [x[0] for x in ret]

            if duplicate_protocol == "replace":
                logger.debug(f"Conflict at {exc.orig.diag.constraint_name}")
                logger.debug(
                    f"Found duplicate entry in {db_name} - "
                    f"{str(exc)}."
                    f"Replacing with a new entry."
                )
                primary_key_val = value_dict[primary_key.name]
                sequence_key_names, sequence_values = new.update_entry(
                    primary_key_val=primary_key_val
                )

        return sequence_key_names, sequence_values

    def execute_query(self, sql_query: str, db_name: str) -> list[Row]:
        """
        Generically execute SQL query

        :param sql_query: SQL query to execute
        :param db_name: db name
        :return: rows from db
        """
        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            logger.debug(f"Query: {sql_query}")

            with conn.execute(sql_query) as cursor:
                query_output = cursor.fetchall()

        return query_output

    @staticmethod
    def create_db(db_name: str):
        """
        Creates a database using credentials, if it does not exist

        :param db_name: DB to create
        :return: None
        """

        engine = get_engine(db_name=db_name)
        if not database_exists(engine.url):
            create_database(engine.url)

        assert database_exists(engine.url)

    @staticmethod
    def has_extension(
        extension_name: str,
        db_name: str,
    ) -> bool:
        """
        Function to create q3c extension and index on table

        :param extension_name: name of extension to check
        :param db_name: name of database to check
        :return: boolean extension exists
        """

        engine = get_engine(db_name=db_name)
        with engine.connect() as conn:
            command = text(
                f"SELECT extname FROM pg_extension WHERE extname='{extension_name}'"
            )
            res = conn.execute(command).all()

        assert len(res) <= 1, "More than one extension found"

        return len(res) == 1
