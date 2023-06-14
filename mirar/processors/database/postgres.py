"""
Module containing postgres util functions
"""
# pylint: disable=not-context-manager
import logging
import os
from glob import glob
from pathlib import Path
from typing import Optional

import numpy as np
import psycopg
from psycopg import errors
from psycopg.rows import Row

from mirar.data import DataBlock
from mirar.errors import ProcessorError
from mirar.processors.database.constraints import DBQueryConstraints
from mirar.processors.database.utils import get_ordered_schema_list

logger = logging.getLogger(__name__)

DB_USER_KEY = "DB_USER"
DB_PASSWORD_KEY = "DB_PWD"

PG_ADMIN_USER_KEY = "PG_ADMIN_USER"
PG_ADMIN_PWD_KEY = "PG_ADMIN_PWD"

DB_USER = os.getenv(DB_USER_KEY)
DB_PASSWORD = os.getenv(DB_PASSWORD_KEY)

ADMIN_USER = os.getenv(PG_ADMIN_USER_KEY, DB_USER)
ADMIN_PASSWORD = os.getenv(PG_ADMIN_PWD_KEY, DB_PASSWORD)

POSTGRES_DUPLICATE_PROTOCOLS = ["fail", "ignore", "replace"]


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

    def create_table(self, schema_path: str | Path, db_name: str):
        """
        Create a database table

        :param schema_path: File to execute
        :param db_name: name of database
        :return: None
        """
        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            with open(schema_path, "r", encoding="utf8") as schema_file:
                conn.execute(schema_file.read())

        logger.info(f"Created table from schema path {schema_path}")

    def create_tables_from_schema(
        self,
        schema_dir: str | Path,
        db_name: str,
    ):
        """
        Creates a db with tables, as described by .sql files in a directory

        :param schema_dir: Directory containing schema files
        :param db_name: name of DB
        :return: None
        """
        schema_files = glob(f"{schema_dir}/*.sql")
        ordered_schema_files = get_ordered_schema_list(schema_files)
        logger.info(f"Creating the following tables - {ordered_schema_files}")
        for schema_file in ordered_schema_files:
            self.create_table(schema_path=schema_file, db_name=db_name)

    def export_to_db(
        self,
        value_dict: dict | DataBlock,
        db_name: str,
        db_table: str,
        duplicate_protocol: str = "fail",
    ) -> tuple[list, list]:
        """
        Export a list of fields in value dict to a batabase table

        :param value_dict: dictionary/DataBlock/other dictonary-like object to export
        :param db_name: name of db to export to
        :param db_table: table of DB to export to
        :param duplicate_protocol: protocol for handling duplicates,
            in "fail"/"ignore"/"replace"
        :return:
        """

        assert duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True

            sql_query = f"""
            SELECT Col.Column_Name from
                INFORMATION_SCHEMA.TABLE_CONSTRAINTS Tab,
                INFORMATION_SCHEMA.CONSTRAINT_COLUMN_USAGE Col
            WHERE
                Col.Constraint_Name = Tab.Constraint_Name
                AND Col.Table_Name = Tab.Table_Name
                AND Constraint_Type = 'PRIMARY KEY'
                AND Col.Table_Name = '{db_table}'
            """
            serial_keys, serial_key_values = [], []
            with conn.execute(sql_query) as cursor:
                primary_key = [x[0] for x in cursor.fetchall()]
                serial_keys = list(self.get_sequence_keys_from_table(db_table, db_name))
                logger.debug(serial_keys)
                colnames = [
                    desc[0]
                    for desc in conn.execute(
                        f"SELECT * FROM {db_table} LIMIT 1"
                    ).description
                    if desc[0] not in serial_keys
                ]

                colnames_str = ""
                for column in colnames:
                    colnames_str += f'"{column}",'
                colnames_str = colnames_str[:-1]
                txt = f"INSERT INTO {db_table} ({colnames_str}) VALUES ("

                for char in ["[", "]", "'"]:
                    txt = txt.replace(char, "")

                for column in colnames:
                    txt += f"'{str(value_dict[column])}', "

                txt = txt + ") "
                txt = txt.replace(", )", ")")

                if len(serial_keys) > 0:
                    txt += "RETURNING "
                    for key in serial_keys:
                        txt += f"{key},"
                    txt += ";"
                    txt = txt.replace(",;", ";")

                logger.debug(txt)
                command = txt

                try:
                    cursor.execute(command)
                    if len(serial_keys) > 0:
                        serial_key_values = cursor.fetchall()[0]
                    else:
                        serial_key_values = []

                except errors.UniqueViolation as exc:
                    primary_key_values = [value_dict[x] for x in primary_key]

                    if duplicate_protocol == "fail":
                        err = (
                            f"Duplicate error, entry with "
                            f"{primary_key}={primary_key_values} "
                            f"already exists in {db_name}."
                        )
                        logger.error(err)
                        raise errors.UniqueViolation from exc

                    if duplicate_protocol == "ignore":
                        logger.debug(
                            f"Found duplicate entry with "
                            f"{primary_key}={primary_key_values} in {db_name}. "
                            f"Ignoring."
                        )
                    elif duplicate_protocol == "replace":
                        logger.debug(
                            f"Updating duplicate entry with "
                            f"{primary_key}={primary_key_values} in {db_name}."
                        )

                        db_constraints = DBQueryConstraints(
                            columns=primary_key,
                            accepted_values=primary_key_values,
                        )

                        update_colnames = []
                        for column in colnames:
                            if column not in primary_key:
                                update_colnames.append(column)

                        if len(serial_keys) == 0:
                            self.modify_db_entry(
                                db_constraints=db_constraints,
                                value_dict=value_dict,
                                db_alter_columns=update_colnames,
                                db_table=db_table,
                                db_name=db_name,
                                return_columns=primary_key,
                            )

                        else:
                            serial_key_values = self.modify_db_entry(
                                db_constraints=db_constraints,
                                value_dict=value_dict,
                                db_alter_columns=update_colnames,
                                db_table=db_table,
                                db_name=db_name,
                                return_columns=primary_key,
                            )

        return serial_keys, serial_key_values

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
        :param db_constraints: constraints to query db
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
                f'"{db_alter_columns[ind]}"=' + f"'{db_alter_values[ind]}'"
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

    def get_sequence_keys_from_table(self, db_table: str, db_name: str) -> np.ndarray:
        """
        Gets sequence keys of db table

        :param db_table: database table to use
        :param db_name: dataname name
        :return: numpy array of keys
        """
        with psycopg.connect(
            f"dbname={db_name} user={self.db_user} password={self.db_password}"
        ) as conn:
            conn.autocommit = True
            sequences = [
                x[0]
                for x in conn.execute(
                    "SELECT c.relname FROM pg_class c WHERE c.relkind = 'S';"
                ).fetchall()
            ]
            seq_tables = np.array([x.split("_")[0] for x in sequences])
            seq_columns = np.array([x.split("_")[1] for x in sequences])
            table_sequence_keys = seq_columns[(seq_tables == db_table)]
        return table_sequence_keys

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
