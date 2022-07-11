import astropy.io.fits
import psycopg
from astropy.io import fits
import os
from glob import glob
from astropy.time import Time
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from psycopg.errors import Error
import logging
from winterdrp.errors import ProcessorError

logger = logging.getLogger(__name__)

default_db_user = os.path.basename(os.environ["HOME"])

schema_dir = os.path.join(os.path.dirname(__file__), "schema")


class DataBaseError(ProcessorError):
    pass


def validate_credentials(
        db_user: str,
        password: str
):
    if db_user is None:
        err = "'db_user' is set as None. Please pass a db_user as an argument, " \
              "or set the environment variable 'PG_DEFAULT_USER'. Using "
        logger.warning(err)
        raise DataBaseError(err)

    if password is None:
        err = "'password' is set as None. Please pass a password as an argument, " \
              "or set the environment variable 'PG_DEFAULT_PWD'."
        logger.error(err)
        raise DataBaseError(err)


def create_db(
        db_name: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD')
):
    validate_credentials(db_user=db_user, password=password)

    with psycopg.connect(f"dbname=postgres user={db_user} password={password}") as conn:
        conn.autocommit = True
        sql = f'''CREATE database {db_name}'''
        conn.execute(sql)
        logger.info(f'Created db {db_name}')


def create_table(
        schema_path: str,
        db_name: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD')
):
    validate_credentials(db_user, password)

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        with open(schema_path, "r") as f:
            conn.execute(f.read())


def create_new_user(
        new_db_user: str,
        new_password: str
):
    default_user = os.environ.get('PG_DEFAULT_USER')
    default_password = os.environ.get('PG_DEFAULT_PWD')

    validate_credentials(new_db_user, new_password)
    validate_credentials(db_user=default_user, password=default_password)

    with psycopg.connect(f"dbname=postgres user={default_user} password={default_password}") as conn:
        conn.autocommit = True
        command = f'''CREATE ROLE {new_db_user} WITH password '{new_password}';'''
        conn.execute(command)


def grant_privileges(
        db_name: str,
        db_user: str
):
    default_user = os.environ.get('PG_DEFAULT_USER')
    default_password = os.environ.get('PG_DEFAULT_PWD')
    validate_credentials(default_user, default_password)

    with psycopg.connect(f"dbname=postgres user={default_user} password={default_password}") as conn:
        conn.autocommit = True
        command = f'''GRANT ALL PRIVILEGES ON DATABASE {db_name} TO {db_user};'''
        conn.execute(command)


def check_if_user_exists(
        user_name: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD')
) -> bool:

    validate_credentials(db_user, password)

    with psycopg.connect(f"dbname=postgres user={db_user} password={password}") as conn:
        conn.autocommit = True
        command = '''SELECT usename FROM pg_user;'''
        data = conn.execute(command).fetchall()
    existing_user_names = [x[0] for x in data]
    logger.debug(f"Found the following users: {existing_user_names}")

    user_exist_bool = user_name in existing_user_names
    logger.info(f"User '{user_name}' {['does not exist', 'already exists'][user_exist_bool]}")
    return user_exist_bool


def check_if_db_exists(
        db_name: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD')
) -> bool:

    validate_credentials(db_user, password)

    with psycopg.connect(f"dbname=postgres user={db_user} password={password}") as conn:
        conn.autocommit = True
        command = '''SELECT datname FROM pg_database;'''
        data = conn.execute(command).fetchall()

    existing_db_names = [x[0] for x in data]
    logger.debug(f"Found the following databases: {existing_db_names}")

    db_exist_bool = db_name in existing_db_names
    logger.info(f"Database '{db_name}' {['does not exist', 'already exists'][db_exist_bool]}")

    return db_exist_bool


def check_if_table_exists(
        db_name: str,
        db_table: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD')
) -> bool:

    validate_credentials(db_user=db_user, password=password)

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        command = '''SELECT datname FROM pg_database;'''
        data = conn.execute(command).fetchall()

    existing_db_names = [x[0] for x in data]
    logger.debug(f"Found the following databases: {existing_db_names}")

    db_exist_bool = db_name in existing_db_names
    logger.info(f"Database table '{db_table}' {['does not exist', 'already exists'][db_exist_bool]}")

    return db_exist_bool


def get_foreign_tables_list(
        schema_files: list[str]
) -> np.ndarray:
    foreign_tables_list = []
    for schema_file in schema_files:
        table_names = []
        with open(schema_file, 'r') as f:
            schema = f.read()
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


def get_ordered_schema_list(
        schema_files: list[str]
) -> list[str]:
    foreign_tables_list = get_foreign_tables_list(schema_files)
    ordered_schema_list = []
    tables_created = []
    schema_table_names = [x.split('/')[-1].split('.sql')[0] for x in schema_files]
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
        schema_dir: str,
        db_name: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD')
):
    schema_files = glob(f'{schema_dir}/*.sql')
    ordered_schema_files = get_ordered_schema_list(schema_files)
    logger.info(f"Creating the following tables - {ordered_schema_files}")
    for schema_file in ordered_schema_files:
        create_table(schema_path=schema_file, db_name=db_name, db_user=db_user, password=password)


def export_to_db(
        value_dict: dict | astropy.io.fits.Header,
        db_name: str,
        db_table: str,
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD'),
) -> tuple[str, list]:
    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
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

        with conn.execute(sql_query) as cursor:

            primary_key = [x[0] for x in cursor.fetchall()]

            colnames = [
                desc[0] for desc in conn.execute(f"SELECT * FROM {db_table} LIMIT 1").description
                if desc[0] not in primary_key
            ]
            logger.debug(primary_key)
            logger.debug(colnames)

            txt = f'INSERT INTO {db_table} ({colnames}) VALUES ('

            for char in ["[", "]", "'"]:
                txt = txt.replace(char, '')

            for c in colnames:
                logger.debug(f"{c}, {value_dict[c.upper()]}")
                txt += f"'{str(value_dict[c.upper()])}', "

            txt = txt + ') '
            txt = txt.replace(', )', ')')
            txt += f"RETURNING "
            for key in primary_key:
                txt += f"{key},"
            txt += ";"
            txt = txt.replace(',;', ';')

            logger.debug(txt)
            command = txt
            cursor.execute(command)
            primary_key_values = cursor.fetchall()[0]

    return primary_key, primary_key_values


def import_from_db(
        db_name: str,
        db_table: str,
        db_query_columns: str | list[str],
        db_accepted_values: str | int | float | list[str | float | int],
        db_output_columns: str | list[str],
        db_user: str = os.environ.get('PG_DEFAULT_USER', default_db_user),
        password: str = os.environ.get('PG_DEFAULT_PWD'),
) -> list[dict]:
    """Query an SQL database with constraints, and return a list of dictionaries.
    One dictionary per entry returned from the query.

    Parameters
    ----------
    db_name: Name of database to query
    db_table: Name of database table to query
    db_query_columns: Name of column to query
    db_accepted_values: Accepted value for query for column
    db_output_columns: Name(s) of columns to return for matched database entries
    db_user: Username for database
    password: password for database

    Returns
    -------
    A list of dictionaries (one per entry)
    """

    if not isinstance(db_query_columns, list):
        db_query_columns = [db_query_columns]

    if not isinstance(db_accepted_values, list):
        db_accepted_values = [db_accepted_values]

    if not isinstance(db_output_columns, list):
        db_output_columns = [db_output_columns]

    assert len(db_query_columns) == len(db_accepted_values)

    all_query_res = []

    constraints = " AND ".join([f"{x} = {db_accepted_values[i]}" for i, x in enumerate(db_query_columns)])

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        sql_query = f"""
        SELECT {', '.join(db_output_columns)} from {db_table}
            WHERE {constraints};
        """

        logger.debug(f"Query: {sql_query}")

        with conn.execute(sql_query) as cursor:
            query_output = cursor.fetchall()

        for entry in query_output:

            assert len(entry) == len(db_output_columns)

            query_res = dict()

            for i, key in enumerate(db_output_columns):
                query_res[key] = entry[i]

            all_query_res.append(query_res)

    return all_query_res
