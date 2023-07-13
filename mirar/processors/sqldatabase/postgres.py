"""
Module containing postgres util functions
"""
# pylint: disable=not-context-manager
import logging
from pathlib import Path
from typing import Optional, Type

import psycopg
from psycopg import errors
from psycopg.rows import Row
from pydantic import ValidationError
from sqlalchemy import inspect
from sqlalchemy.exc import IntegrityError

from mirar.data import DataBlock
from mirar.errors import ProcessorError
from mirar.processors.database.constraints import DBQueryConstraints
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.processors.sqldatabase.postgres_utils import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_PASSWORD_KEY,
    DB_USER,
    DB_USER_KEY,
    PG_ADMIN_PWD_KEY,
    PG_ADMIN_USER_KEY,
    POSTGRES_DUPLICATE_PROTOCOLS,
)

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

        # TODO check user exists

    def run_sql_command_from_file(self, file_path: str | Path, db_name: str):
        """
        Execute SQL command from file

        :param file_path: File to execute
        :param db_name: name of database
        :return: False
        """
        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            with open(file_path, "r", encoding="utf8") as sql_file:
                conn.execute(sql_file.read())

            logger.info(f"Executed sql commands from file {file_path}")

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

    def import_from_db(
        self,
        db_name: str,
        db_table: str,
        db_output_columns: str | list[str],
        output_alias_map: Optional[str | list[str]] = None,
        max_num_results: Optional[int] = None,
        db_constraints: Optional[DBQueryConstraints] = None,
    ) -> list[dict]:
        """Query an SQL database with constraints, and return a list of dictionaries.
        One dictionary per entry returned from the query.

        #TODO check admin

        Parameters
        ----------
        db_name: Name of database to query
        db_table: Name of database table to query
        db_output_columns: Name(s) of columns to return for matched database entries
        output_alias_map: Alias to assign for each output column
        max_num_results: Maximum number of results to return

        Returns
        -------
        A list of dictionaries (one per entry)
        """

        if not isinstance(db_output_columns, list):
            db_output_columns = [db_output_columns]

        if output_alias_map is None:
            output_alias_map = db_output_columns

        if not isinstance(output_alias_map, list):
            output_alias_map = [output_alias_map]

        assert len(output_alias_map) == len(db_output_columns)

        all_query_res = []

        if db_constraints is not None:
            constraints = db_constraints.parse_constraints()
        else:
            constraints = ""

        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            sql_query = f"""
            SELECT {', '.join(db_output_columns)} from {db_table}
                WHERE {constraints}
            """

            if max_num_results is not None:
                sql_query += f" LIMIT {max_num_results}"

            sql_query += ";"

            logger.debug(f"Query: {sql_query}")

            with conn.execute(sql_query) as cursor:
                query_output = cursor.fetchall()

            for entry in query_output:
                assert len(entry) == len(db_output_columns)

                query_res = {}

                for i, key in enumerate(output_alias_map):
                    query_res[key] = entry[i]

                all_query_res.append(query_res)

        return all_query_res

    def crossmatch_with_database(
        self,
        db_name: str,
        db_table: str,
        db_output_columns: str | list[str],
        ra: float,
        dec: float,
        crossmatch_radius_arcsec: float,
        output_alias_map: Optional[dict] = None,
        ra_field_name: str = "ra",
        dec_field_name: str = "dec",
        query_distance_bool: bool = False,
        q3c_bool: bool = False,
        query_constraints: Optional[DBQueryConstraints] = None,
        order_field_name: Optional[str] = None,
        num_limit: Optional[int] = None,
    ) -> list[dict]:
        """
        Crossmatch a given spatial position (ra/dec) with sources in a database,
        and returns a list of matches

        #TODO: check admin

        :param db_name: name of db to query
        :param db_table: name of db table
        :param db_output_columns: columns to return
        :param output_alias_map: mapping for renaming columns
        :param ra: RA
        :param dec: dec
        :param crossmatch_radius_arcsec: radius for crossmatch
        :param ra_field_name: name of ra column in database
        :param dec_field_name: name of dec column in database
        :param query_distance_bool: boolean where to return crossmatch distance
        :param q3c_bool: boolean whether to use q3c_bool
        :param order_field_name: field to order result by
        :param num_limit: limit on sql query
        :return: list of query result dictionaries
        """

        if output_alias_map is None:
            output_alias_map = {}
            for col in db_output_columns:
                output_alias_map[col] = col

        crossmatch_radius_deg = crossmatch_radius_arcsec / 3600.0

        if q3c_bool:
            constraints = (
                f"q3c_radial_query({ra_field_name},{dec_field_name},"
                f"{ra},{dec},{crossmatch_radius_deg}) "
            )
        else:
            ra_min = ra - crossmatch_radius_deg
            ra_max = ra + crossmatch_radius_deg
            dec_min = dec - crossmatch_radius_deg
            dec_max = dec + crossmatch_radius_deg
            constraints = (
                f" {ra_field_name} between {ra_min} and {ra_max} AND "
                f"{dec_field_name} between {dec_min} and {dec_max} "
            )

        if query_constraints is not None:
            constraints += f"""AND {query_constraints.parse_constraints()}"""

        select = f""" {'"' + '","'.join(db_output_columns) + '"'}"""
        if query_distance_bool:
            if q3c_bool:
                select = (
                    f"q3c_dist({ra_field_name},{dec_field_name},{ra},{dec}) AS xdist,"
                    + select
                )
            else:
                select = f"""{ra_field_name} - ra AS xdist,""" + select

        query = f"""SELECT {select} FROM {db_table} WHERE {constraints}"""

        if order_field_name is not None:
            query += f""" ORDER BY {order_field_name}"""
        if num_limit is not None:
            query += f""" LIMIT {num_limit}"""

        query += ";"

        query_output = self.execute_query(query, db_name)
        all_query_res = []

        for entry in query_output:
            if not query_distance_bool:
                assert len(entry) == len(db_output_columns)
            else:
                assert len(entry) == len(db_output_columns) + 1
            query_res = {}
            for i, key in enumerate(output_alias_map):
                query_res[key] = entry[i]
                if query_distance_bool:
                    query_res["dist"] = entry["xdist"]
            all_query_res.append(query_res)
        return all_query_res

    def check_if_exists(
        self, check_command: str, check_value: str, db_name: str = "postgres"
    ) -> bool:
        """
        Check if a user account exists

        :param check_command if a user/database/table exists
        :param check_value: username to check
        :param db_name: name of database to query
        :return: boolean
        """
        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            data = conn.execute(check_command).fetchall()
        existing_user_names = [x[0] for x in data]
        logger.debug(f"Found the following values: {existing_user_names}")

        return check_value in existing_user_names

    def check_if_db_exists(self, db_name: str) -> bool:
        """
        Check if a user account exists

        :param db_name: database to check
        :return: boolean
        """

        check_command = """SELECT datname FROM pg_database;"""

        db_exist_bool = self.check_if_exists(
            check_command=check_command,
            check_value=db_name,
            db_name="postgres",
        )

        logger.debug(f"Database '{db_name}' does {['not ', ''][db_exist_bool]} exist")

        return db_exist_bool

    def create_db(self, db_name: str):
        """
        Creates a database using credentials

        :param db_name: DB to create
        :return: None
        """

        with psycopg.connect(
            f"dbname=postgres user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            sql = f"""CREATE database {db_name}"""
            conn.execute(sql)
            logger.info(f"Created db {db_name}")

    def check_if_table_exists(self, db_name: str, db_table: str) -> bool:
        """
        Check if a db table account exists

        :param db_name: database to check
        :param db_table: table to check
        :return: boolean
        """

        check_command = (
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='public';"
        )

        table_exist_bool = self.check_if_exists(
            check_command=check_command,
            check_value=db_table,
            db_name=db_name,
        )

        logger.debug(f"Table '{db_table}' does {['not ', ''][table_exist_bool]} exist")

        return table_exist_bool

    def modify_db_entry(
        self,
        db_name: str,
        db_table: str,
        db_constraints: DBQueryConstraints,
        value_dict: dict | DataBlock,
        db_alter_columns: str | list[str],
        return_columns: Optional[str | list[str]] = None,
    ) -> list[Row]:
        """
        Modify a db entry

        :param db_name: name of db
        :param db_table: Name of table
        :param value_dict: dict-like object to provide updated values
        :param db_alter_columns: columns to alter in db
        :param return_columns: columns to return
        :return: db query (return columns)
        """

        if not isinstance(db_alter_columns, list):
            db_alter_columns = [db_alter_columns]

        if return_columns is None:
            return_columns = db_alter_columns
        if not isinstance(return_columns, list):
            return_columns = [return_columns]

        constraints = db_constraints.parse_constraints()

        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True

            db_alter_values = [str(value_dict[c]) for c in db_alter_columns]

            alter_values_txt = [
                f"{db_alter_columns[ind]}='{db_alter_values[ind]}'"
                for ind in range(len(db_alter_columns))
            ]

            sql_query = (
                f"UPDATE {db_table} SET {', '.join(alter_values_txt)} "
                f"WHERE {constraints}"
            )

            if len(return_columns) > 0:
                logger.debug(return_columns)
                sql_query += f""" RETURNING {', '.join(return_columns)}"""
            sql_query += ";"
            query_output = self.execute_query(sql_query, db_name)

        return query_output


class PostgresAdmin(PostgresUser):
    """
    An Admin postgres user, with additional functionality for creatying new users
    """

    user_env_varaiable = PG_ADMIN_USER_KEY
    pass_env_variable = PG_ADMIN_PWD_KEY

    def __init__(self, db_user: str = ADMIN_USER, db_password: str = ADMIN_PASSWORD):
        super().__init__(db_user=db_user, db_password=db_password)

    def create_new_user(self, new_db_user: str, new_password: str):
        """
        Create a new postgres user

        :param new_db_user: new username
        :param new_password: new user password
        :return: None
        """

        with psycopg.connect(
            f"dbname=postgres user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            command = f"CREATE ROLE {new_db_user} WITH password '{new_password}' LOGIN;"
            conn.execute(command)

    def grant_privileges(self, db_name: str, db_user: str, schema: str = "public"):
        """
        Grant privilege to user on database

        :param db_name: name of database
        :param db_user: username to grant privileges for db_user
        :return: None
        """
        with psycopg.connect(
            f"dbname=postgres user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            command = f"""GRANT ALL PRIVILEGES ON DATABASE {db_name} TO {db_user};"""
            conn.execute(command)
            command = f"""ALTER DATABASE {db_name} OWNER TO {db_user}"""
            conn.execute(command)

        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            command = f"""GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA {schema}
            TO {db_user};"""
            conn.execute(command)
            command = f"""GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA {schema}
                        TO {db_user};"""
            conn.execute(command)

    def check_if_user_exists(self, user_name: str) -> bool:
        """
        Check if a user account exists

        :param user_name: username to check
        :return: boolean
        """
        check_command = """SELECT usename FROM pg_user;"""

        user_exist_bool = self.check_if_exists(
            check_command=check_command,
            check_value=user_name,
        )

        logger.debug(f"User '{user_name}' does {['not ', ''][user_exist_bool]} exist")

        return user_exist_bool
