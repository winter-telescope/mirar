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

from winterdrp.data import DataBlock
from winterdrp.errors import ProcessorError
from winterdrp.processors.database.constraints import DBQueryConstraints

logger = logging.getLogger(__name__)

PG_ADMIN_USER_KEY = "PG_ADMIN_USER"
PG_ADMIN_PWD_KEY = "PG_ADMIN_PWD"

POSTGRES_DUPLICATE_PROTOCOLS = ["fail", "ignore", "replace"]


class DataBaseError(ProcessorError):
    """Error relating to postgres interactions"""


def validate_credentials(db_user: str, password: str, admin: bool = False):
    """
    Checks that user credentials exist

    :param db_user: Username
    :param password: Password
    :param admin: boolean whether the account is an admin one
    :return: None
    """

    if db_user is None:
        if admin:
            user = "admin_db_user"
            env_user_var = PG_ADMIN_USER_KEY
        else:
            user = "db_user"
            env_user_var = "DB_USER"
        err = (
            f"'{user}' is set as None. Please pass a db_user as an argument, "
            f"or set the environment variable '{env_user_var}'. Using "
        )
        logger.warning(err)
        raise DataBaseError(err)

    if password is None:
        if admin:
            pwd = "db_admin_password"
            env_pwd_var = PG_ADMIN_PWD_KEY
        else:
            pwd = "password"
            env_pwd_var = "DB_PWD"
        err = (
            f"'{pwd}' is set as None. Please pass a password as an argument, "
            f"or set the environment variable '{env_pwd_var}'."
        )
        logger.error(err)
        raise DataBaseError(err)


def create_db(db_name: str):
    """
    Creates a database using credentials

    :param db_name: DB to create
    :return: None
    """
    admin_user = os.environ.get(PG_ADMIN_USER_KEY)
    admin_password = os.environ.get(PG_ADMIN_PWD_KEY)

    validate_credentials(db_user=admin_user, password=admin_password)

    with psycopg.connect(
        f"dbname=postgres user={admin_user} password={admin_password}"
    ) as conn:
        conn.autocommit = True
        sql = f"""CREATE database {db_name}"""
        conn.execute(sql)
        logger.info(f"Created db {db_name}")


def run_sql_command_from_file(
    file_path: str | Path,
    db_name: str,
    db_user: str,
    password: str,
    admin: bool = False,
):
    """
    Execute SQL command from file

    :param file_path: File to execute
    :param db_name: name of database
    :param db_user: Postgres db user
    :param password: Postgress password
    :param admin: Whether to use an admin user
    :return: False
    """
    validate_credentials(db_name, db_user, admin)
    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={password}"
    ) as conn:
        with open(file_path, "r", encoding="utf8") as sql_file:
            conn.execute(sql_file.read())

        logger.info(f"Executed sql commands from file {file_path}")


def create_table(schema_path: str | Path, db_name: str, db_user: str, password: str):
    """
    Create a database table

    :param schema_path: File to execute
    :param db_name: name of database
    :param db_user: Postgres db user
    :param password: Postgress password
    :return: False
    """
    validate_credentials(db_user, password)

    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={password}"
    ) as conn:
        conn.autocommit = True
        with open(schema_path, "r", encoding="utf8") as schema_file:
            conn.execute(schema_file.read())

    logger.info(f"Created table from schema path {schema_path}")


def create_new_user(new_db_user: str, new_password: str):
    """
    Create a new postgres user

    :param new_db_user: new username
    :param new_password: new user password
    :return: None
    """
    admin_user = os.environ.get(PG_ADMIN_USER_KEY)
    admin_password = os.environ.get(PG_ADMIN_PWD_KEY)

    validate_credentials(new_db_user, new_password)
    validate_credentials(db_user=admin_user, password=admin_password, admin=True)

    with psycopg.connect(
        f"dbname=postgres user={admin_user} password={admin_password}"
    ) as conn:
        conn.autocommit = True
        command = f"""CREATE ROLE {new_db_user} WITH password '{new_password}' LOGIN;"""
        conn.execute(command)


def grant_privileges(db_name: str, db_user: str):
    """
    Grant privilege to user on database

    :param db_name: name of database
    :param db_user: username
    :return: None
    """
    admin_user = os.environ.get(PG_ADMIN_USER_KEY)
    admin_password = os.environ.get(PG_ADMIN_PWD_KEY)
    validate_credentials(admin_user, admin_password, admin=True)

    with psycopg.connect(
        f"dbname=postgres user={admin_user} password={admin_password}"
    ) as conn:
        conn.autocommit = True
        command = f"""GRANT ALL PRIVILEGES ON DATABASE {db_name} TO {db_user};"""
        conn.execute(command)


def check_if_exists(
    check_command: str,
    check_value: str,
    db_name: str = "postgres",
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    password: str = os.environ.get(PG_ADMIN_PWD_KEY),
) -> bool:
    """
    Check if a user account exists

    :param check_command if a user/database/table exists
    :param check_value: username to check
    :param db_name: name of database to query
    :param db_user: username to use for check
    :param password: password to use for check
    :return: boolean
    """
    validate_credentials(db_user, password)

    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={password}"
    ) as conn:
        conn.autocommit = True
        data = conn.execute(check_command).fetchall()
    existing_user_names = [x[0] for x in data]
    logger.debug(f"Found the following values: {existing_user_names}")

    return check_value in existing_user_names


def check_if_user_exists(
    user_name: str,
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    password: str = os.environ.get(PG_ADMIN_PWD_KEY),
) -> bool:
    """
    Check if a user account exists

    :param user_name: username to check
    :param db_user: username to use for check
    :param password: password to use for check
    :return: boolean
    """
    check_command = """SELECT usename FROM pg_user;"""

    user_exist_bool = check_if_exists(
        check_command=check_command,
        check_value=user_name,
        db_user=db_user,
        password=password,
    )

    logger.debug(
        f"User '{user_name}' {['does not exist', 'already exists'][user_exist_bool]}"
    )
    return user_exist_bool


def check_if_db_exists(
    db_name: str,
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    password: str = os.environ.get(PG_ADMIN_PWD_KEY),
) -> bool:
    """
    Check if a user account exists

    :param db_name: database to check
    :param db_user: username to use for check
    :param password: password to use for check
    :return: boolean
    """

    check_command = """SELECT datname FROM pg_database;"""

    db_exist_bool = check_if_exists(
        check_command=check_command,
        check_value=db_name,
        db_user=db_user,
        password=password,
        db_name="postgres",
    )

    logger.debug(
        f"Database '{db_name}' {['does not exist', 'already exists'][db_exist_bool]}"
    )

    return db_exist_bool


def check_if_table_exists(
    db_name: str, db_table: str, db_user: str, password: str
) -> bool:
    """
    Check if a db table account exists

    :param db_name: database to check
    :param db_table: table to check
    :param db_user: username to use for check
    :param password: password to use for check
    :return: boolean
    """

    check_command = (
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema='public';"
    )

    table_exist_bool = check_if_exists(
        check_command=check_command,
        check_value=db_table,
        db_name=db_name,
        db_user=db_user,
        password=password,
    )

    logger.debug(
        f"Database table '{db_table}' "
        f"{['does not exist', 'already exists'][table_exist_bool]}"
    )

    return table_exist_bool


def get_foreign_tables_list(schema_files: list[str]) -> np.ndarray:
    """
    Returns a list of foreign tables

    :param schema_files: List of schema files to read
    :return: Returns list of foreign keys in schema
    """
    foreign_tables_list = []
    for schema_file_path in schema_files:
        table_names = []
        with open(schema_file_path, "r", encoding="utf8") as schema_file:
            schema = schema_file.read()
        if "FOREIGN KEY" not in schema:
            pass
        else:
            schema = schema.replace("\n", "")
            schema = schema.replace("\t", "")
            schema_split = np.array(schema.split(","))
            fk_rows = np.array(["FOREIGN KEY" in x for x in schema_split])
            for row in schema_split[fk_rows]:
                words = np.array(row.split(" "))
                refmask = np.array(["REFERENCES" in x for x in words])
                idx = np.where(refmask)[0][0] + 1
                tablename = words[idx].split("(")[0]
                table_names.append(tablename)
        foreign_tables_list.append(np.array(table_names))
    return np.array(foreign_tables_list)


def get_ordered_schema_list(schema_files: list[str]) -> list[str]:
    """
    Returns an ordered list of schema, ensuring connected database tables
    are created in the right order

    :param schema_files: List of schema files to read
    :return: Returns list of foreign keys in schema
    """
    foreign_tables_list = get_foreign_tables_list(schema_files)
    ordered_schema_list = []
    tables_created = []
    schema_table_names = [x.split("/")[-1].split(".sql")[0] for x in schema_files]
    while len(tables_created) < len(schema_files):
        for ind, schema_file in enumerate(schema_files):
            table_name = schema_table_names[ind]
            if table_name in tables_created:
                pass
            else:
                foreign_tables = foreign_tables_list[ind]
                if len(foreign_tables) == 0:
                    ordered_schema_list.append(schema_file)
                    tables_created.append(table_name)
                else:
                    if np.all(np.isin(foreign_tables, tables_created)):
                        ordered_schema_list.append(schema_file)
                        tables_created.append(table_name)

    return ordered_schema_list


def create_tables_from_schema(
    schema_dir: str | Path,
    db_name: str,
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    password: str = os.environ.get(PG_ADMIN_PWD_KEY),
):
    """
    Creates a db with tables, as described by .sql files in a directory

    :param schema_dir: Directory containing schema files
    :param db_name: name of DB
    :param db_user: db user
    :param password: db password
    :return: None
    """
    schema_files = glob(f"{schema_dir}/*.sql")
    ordered_schema_files = get_ordered_schema_list(schema_files)
    logger.info(f"Creating the following tables - {ordered_schema_files}")
    for schema_file in ordered_schema_files:
        create_table(
            schema_path=schema_file, db_name=db_name, db_user=db_user, password=password
        )


def export_to_db(
    value_dict: dict | DataBlock,
    db_name: str,
    db_table: str,
    db_user: str,
    password: str,
    duplicate_protocol: str = "fail",
) -> tuple[list, list]:
    """
    Export a list of fields in value dict to a batabase table

    :param value_dict: dictionary/DataBlock/other dictonary-like object to export
    :param db_name: name of db to export to
    :param db_table: table of DB to export to
    :param db_user: db user
    :param password: password
    :param duplicate_protocol: protocol for handling duplicates,
        in "fail"/"ignore"/"replace"
    :return:
    """

    assert duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={password}"
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
            serial_keys = list(
                get_sequence_keys_from_table(db_table, db_name, db_user, password)
            )
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
                        f"{primary_key}={primary_key_values} in {db_name}. Ignoring."
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

                    serial_key_values = modify_db_entry(
                        db_constraints=db_constraints,
                        value_dict=value_dict,
                        db_alter_columns=update_colnames,
                        db_table=db_table,
                        db_name=db_name,
                        db_user=db_user,
                        password=password,
                        return_columns=serial_keys,
                    )

    return serial_keys, serial_key_values


def import_from_db(
    db_name: str,
    db_table: str,
    db_output_columns: str | list[str],
    output_alias_map: Optional[str | list[str]] = None,
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    password: str = os.environ.get(PG_ADMIN_PWD_KEY),
    max_num_results: Optional[int] = None,
    db_constraints: Optional[DBQueryConstraints] = None,
) -> list[dict]:
    """Query an SQL database with constraints, and return a list of dictionaries.
    One dictionary per entry returned from the query.

    Parameters
    ----------
    db_name: Name of database to query
    db_table: Name of database table to query
    db_output_columns: Name(s) of columns to return for matched database entries
    output_alias_map: Alias to assign for each output column
    db_user: Username for database
    password: password for database
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
        f"dbname={db_name} user={db_user} password={password}"
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


def execute_query(
    sql_query: str, db_name: str, db_user: str, password: str
) -> list[Row]:
    """
    Generically execute SQL query

    :param sql_query: SQL query to execute
    :param db_name: db name
    :param db_user: db user
    :param password: db password
    :return: rows from db
    """
    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={password}"
    ) as conn:
        conn.autocommit = True
        logger.debug(f"Query: {sql_query}")

        with conn.execute(sql_query) as cursor:
            query_output = cursor.fetchall()

    return query_output


def crossmatch_with_database(
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
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    db_password: str = os.environ.get(PG_ADMIN_PWD_KEY),
) -> list[dict]:
    """
    Crossmatch a given spatial position (ra/dec) with sources in a database,
    and returns a list of matches

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
    :param db_user: db user
    :param db_password: db password
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
                f"""q3c_dist({ra_field_name},{dec_field_name},{ra},{dec}) AS xdist,"""
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

    query_output = execute_query(query, db_name, db_user, db_password)
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


def get_sequence_keys_from_table(
    db_table: str, db_name: str, db_user: str, password: str
) -> np.ndarray:
    """
    Gets sequence keys of db table

    :param db_table: database table to use
    :param db_name: dataname name
    :param db_user: db user
    :param password: db password
    :return: numpy array of keys
    """
    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={password}"
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


def modify_db_entry(
    db_name: str,
    db_table: str,
    db_constraints: DBQueryConstraints,
    value_dict: dict | DataBlock,
    db_alter_columns: str | list[str],
    return_columns: Optional[str | list[str]] = None,
    db_user: str = os.environ.get(PG_ADMIN_USER_KEY),
    password: str = os.environ.get(PG_ADMIN_PWD_KEY),
) -> list[Row]:
    """
    Modify a db entry

    :param db_name: name of db
    :param db_table: Name of table
    :param value_dict: dict-like object to provide updated values
    :param db_alter_columns: columns to alter in db
    :param return_columns: columns to return
    :param db_user: db user
    :param password: db password
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
        f"dbname={db_name} user={db_user} password={password}"
    ) as conn:
        conn.autocommit = True

        db_alter_values = [str(value_dict[c]) for c in db_alter_columns]

        alter_values_txt = [
            f"{db_alter_columns[ind]}='{db_alter_values[ind]}'"
            for ind in range(len(db_alter_columns))
        ]

        sql_query = f"""
                    UPDATE {db_table} SET {', '.join(alter_values_txt)} WHERE {constraints}
                    """
        if len(return_columns) > 0:
            logger.debug(return_columns)
            sql_query += f""" RETURNING {', '.join(return_columns)}"""
        sql_query += ";"
        query_output = execute_query(sql_query, db_name, db_user, password)

    return query_output


def get_column_names_from_schema(schema_file_path: str | Path) -> list[str]:
    """
    Get column names from a schema file

    :param schema_file_path: file to read
    :return: list of columns
    """
    with open(schema_file_path, "r", encoding="utf8") as schema_file:
        dat = schema_file.read()
    dat = dat.split(");")[0]
    dat = dat.split("\n")[1:-1]
    pkstrip = [x.strip(",").split("PRIMARY KEY")[0].strip() for x in dat]
    fkstrip = [x.strip(",").split("FOREIGN KEY")[0].strip() for x in pkstrip]
    colnames = [x.split(" ")[0].strip('"') for x in fkstrip]
    return colnames
