import astropy.io.fits
import psycopg
import os
from glob import glob
import numpy as np
import logging
from winterdrp.errors import ProcessorError

logger = logging.getLogger(__name__)

schema_dir = os.path.join(os.path.dirname(__file__), "schema")

pg_admin_user_key = 'PG_ADMIN_USER'
pg_admin_pwd_key = 'PG_ADMIN_PWD'

class DataBaseError(ProcessorError):
    pass


def validate_credentials(
        db_user: str,
        password: str,
        admin=False
):

    if db_user is None:
        user = 'db_user'
        env_user_var = 'DB_USER'
        if admin:
            user = 'admin_db_user'
            env_user_var = pg_admin_user_key
        err = f"'{user}' is set as None. Please pass a db_user as an argument, " \
              f"or set the environment variable '{env_user_var}'. Using "
        logger.warning(err)
        raise DataBaseError(err)

    if password is None:
        pwd = 'password'
        env_pwd_var = 'DB_PWD'
        if admin:
            pwd = 'db_admin_password'
            env_pwd_var = pg_admin_pwd_key
        err = f"'{pwd}' is set as None. Please pass a password as an argument, " \
              f"or set the environment variable '{env_pwd_var}'."
        logger.error(err)
        raise DataBaseError(err)


def create_db(
        db_name: str,
):
    admin_user = os.environ.get(pg_admin_user_key)
    admin_password = os.environ.get(pg_admin_pwd_key)
    validate_credentials(db_user=admin_user, password=admin_password)

    with psycopg.connect(f"dbname=postgres user={admin_user} password={admin_password}") as conn:
        conn.autocommit = True
        sql = f'''CREATE database {db_name}'''
        conn.execute(sql)
        logger.info(f'Created db {db_name}')


def create_table(
        schema_path: str,
        db_name: str,
        db_user: str,
        password: str
):
    validate_credentials(db_user, password)

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        with open(schema_path, "r") as f:
            conn.execute(f.read())

    logger.info(f'Created table from schema path {schema_path}')


def create_new_user(
        new_db_user: str,
        new_password: str
):
    admin_user = os.environ.get(pg_admin_user_key)
    admin_password = os.environ.get(pg_admin_pwd_key)

    validate_credentials(new_db_user, new_password)
    validate_credentials(db_user=admin_user, password=admin_password,admin=True)

    with psycopg.connect(f"dbname=postgres user={admin_user} password={admin_password}") as conn:
        conn.autocommit = True
        command = f'''CREATE ROLE {new_db_user} WITH password '{new_password}' LOGIN;'''
        conn.execute(command)


def grant_privileges(
        db_name: str,
        db_user: str
):
    admin_user = os.environ.get(pg_admin_user_key)
    admin_password = os.environ.get(pg_admin_pwd_key)
    validate_credentials(admin_user, admin_password, admin=True)

    with psycopg.connect(f"dbname=postgres user={admin_user} password={admin_password}") as conn:
        conn.autocommit = True
        command = f'''GRANT ALL PRIVILEGES ON DATABASE {db_name} TO {db_user};'''
        conn.execute(command)


def check_if_user_exists(
        user_name: str,
        db_user: str = os.environ.get(pg_admin_user_key),
        password: str = os.environ.get(pg_admin_pwd_key)
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
        db_user: str = os.environ.get(pg_admin_user_key),
        password: str = os.environ.get(pg_admin_pwd_key)
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
        db_user: str ,
        password: str
) -> bool:
    validate_credentials(db_user=db_user, password=password)

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        command = '''SELECT table_name FROM information_schema.tables WHERE table_schema='public';'''
        data = conn.execute(command).fetchall()

    existing_table_names = [x[0] for x in data]
    logger.debug(f"Found the following tables: {existing_table_names}")

    table_exist_bool = db_table in existing_table_names
    logger.info(f"Database table '{db_table}' {['does not exist', 'already exists'][table_exist_bool]}")

    return table_exist_bool


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
        db_user: str = os.environ.get(pg_admin_user_key),
        password: str = os.environ.get(pg_admin_pwd_key)
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
        db_user: str = os.environ.get(pg_admin_user_key),
        password: str = os.environ.get(pg_admin_pwd_key),
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
            sequences = [x[0] for x in conn.execute(f"SELECT c.relname FROM pg_class c WHERE c.relkind = 'S';").fetchall()]
            seq_tables = np.array([x.split('_')[0] for x in sequences])
            seq_columns = np.array([x.split('_')[1] for x in sequences])
            serial_keys = seq_columns[(seq_tables==db_table)]
            logger.debug(serial_keys)
            colnames = [
                desc[0] for desc in conn.execute(f"SELECT * FROM {db_table} LIMIT 1").description
                if desc[0] not in serial_keys
            ]
            logger.debug(primary_key)
            logger.debug(colnames)

            colnames_str = ''
            for x in colnames:
                colnames_str += f'"{x}",'
            colnames_str = colnames_str[:-1]
            txt = f'INSERT INTO {db_table} ({colnames_str}) VALUES ('

            for char in ["[", "]", "'"]:
                txt = txt.replace(char, '')

            for c in colnames:
                logger.debug(f"{c}, {value_dict[c]}")
                txt += f"'{str(value_dict[c])}', "

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


def parse_constraints(db_query_columns,
                      db_comparison_types,
                      db_accepted_values
                      ):
    assert len(db_comparison_types) == len(db_accepted_values)
    assert np.all(np.isin(np.unique(db_comparison_types), ['=', '<', '>', 'between']))
    constraints = ""
    for i, x in enumerate(db_query_columns):
        if db_comparison_types[i] == 'between':
            assert len(db_accepted_values[i]) == 2
            constraints += f"{x} between {db_accepted_values[i][0]} and {db_accepted_values[i][1]} AND "
        else:
            constraints += f"{x} {db_comparison_types[i]} {db_accepted_values[i]} AND "

        constraints = constraints[:-4]  # strip the last AND

    return constraints


def parse_select():
    pass


def import_from_db(
        db_name: str,
        db_table: str,
        db_query_columns: str | list[str],
        db_accepted_values: str | int | float | list[str | float | int | list],
        db_output_columns: str | list[str],
        output_alias_map: str | list[str],
        db_user: str = os.environ.get(pg_admin_user_key),
        password: str = os.environ.get(pg_admin_pwd_key),
        max_num_results: int = None,
        db_comparison_types: list[str] = None
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
    output_alias_map: Alias to assign for each output column
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

    if output_alias_map is None:
        output_alias_map = db_output_columns

    if not isinstance(output_alias_map, list):
        output_alias_map = [output_alias_map]

    assert len(output_alias_map) == len(db_output_columns)

    all_query_res = []

    if db_comparison_types is None:
        db_comparison_types = ['='] * len(db_accepted_values)
    assert len(db_comparison_types) == len(db_accepted_values)
    assert np.isin(np.all(np.unique(db_comparison_types), ['=', '<', '>', 'between']))

    constraints = parse_constraints(db_query_columns,
                                    db_comparison_types,
                                    db_accepted_values)
    # constraints = " AND ".join([f"{x} {db_comparison_types[i]} {db_accepted_values[i]}" for i, x in enumerate(
    # db_query_columns)])

    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        sql_query = f"""
        SELECT {', '.join(db_output_columns)} from {db_table}
            WHERE {constraints}
        """

        if max_num_results is not None:
            sql_query += f" LIMIT {max_num_results}"

        sql_query += f";"

        logger.debug(f"Query: {sql_query}")

        with conn.execute(sql_query) as cursor:
            query_output = cursor.fetchall()

        for entry in query_output:

            assert len(entry) == len(db_output_columns)

            query_res = dict()

            for i, key in enumerate(output_alias_map):
                query_res[key] = entry[i]

            all_query_res.append(query_res)

    return all_query_res


def execute_query(sql_query, db_name, db_user, password):
    with psycopg.connect(f"dbname={db_name} user={db_user} password={password}") as conn:
        conn.autocommit = True
        logger.debug(f"Query: {sql_query}")

        with conn.execute(sql_query) as cursor:
            query_output = cursor.fetchall()

        return query_output


def xmatch_import_db(db_name: str,
                     db_table: str,
                     db_query_columns: str | list[str],
                     db_accepted_values: str | int | float | list[str | float | int],
                     db_output_columns: str | list[str],
                     output_alias_map: str | list[str],
                     ra: float,
                     dec: float,
                     xmatch_radius_arcsec: float,
                     ra_field_name: str = 'ra',
                     dec_field_name: str = 'dec',
                     query_dist=False,
                     q3c=False,
                     db_comparison_types: list[str] = None,
                     order_field_name: str = None,
                     order_ascending: bool = True,
                     num_limit: int = None,
                     db_user: str = os.environ.get(pg_admin_user_key),
                     db_password: str = os.environ.get(pg_admin_pwd_key),
                     ) -> list[dict]:

    if output_alias_map is None:
        output_alias_map = {}
        for col in db_output_columns:
            output_alias_map[col] = col

    xmatch_radius_deg = xmatch_radius_arcsec / 3600.0

    if q3c:
        constraints = f"""q3c_radial_query({ra_field_name},{dec_field_name},{ra},{dec},{xmatch_radius_deg}) """ \

    else:
        ra_min = ra - xmatch_radius_deg
        ra_max = ra + xmatch_radius_deg
        dec_min = dec - xmatch_radius_deg
        dec_max = dec + xmatch_radius_deg
        constraints = f""" {ra_field_name} between {ra_min} and {ra_max} AND {dec_field_name} between {dec_min} and {dec_max} """

    parsed_constraints = parse_constraints(db_query_columns,
                                           db_comparison_types,
                                           db_accepted_values)
    if len(parsed_constraints) > 0:
        constraints += f"""AND {constraints}"""

    select = f""" {'"' + '","'.join(db_output_columns) + '"'}"""
    if query_dist:
        if q3c:
            select = f"""q3c_dist({ra_field_name},{dec_field_name},{ra},{dec}) AS xdist,""" + select
        else:
            select = f"""{ra_field_name} - ra AS xdist,""" + select

    query = f"""SELECT {select} FROM {db_table} WHERE {constraints}"""
    order_seq = ["asc", "desc"][np.sum(order_ascending)]
    if order_field_name is not None:
        query += f""" ORDER BY {order_field_name}"""
    if num_limit is not None:
        query += f""" LIMIT {num_limit}"""

    query += ";"

    query_output = execute_query(query, db_name, db_user, db_password)
    all_query_res = []

    for entry in query_output:
        if not query_dist:
            assert len(entry) == len(db_output_columns)
        else:
            assert len(entry) == len(db_output_columns) + 1
        query_res = dict()
        for i, key in enumerate(output_alias_map):
            query_res[key] = entry[i]
            if query_dist:
                query_res['dist'] = entry['xdist']
        all_query_res.append(query_res)
    return all_query_res


def get_colnames_from_schema(schema_file):
    with open(schema_file,'r') as f:
        dat = f.read()
    dat = dat.split(');')[0]
    dat = dat.split('\n')[1:-1]
    pkstrip = [x.strip(',').split('PRIMARY KEY')[0].strip() for x in dat]
    fkstrip = [x.strip(',').split('FOREIGN KEY')[0].strip() for x in pkstrip]
    colnames = [x.split(' ')[0].strip('"') for x in fkstrip]
    return colnames