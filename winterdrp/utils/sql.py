"""
Util functions for database interactions
"""
from sqlalchemy import DDL, Engine, create_engine

from winterdrp.processors.database.postgres import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_USER,
)


def get_engine(
    db_user: str = DB_USER,
    db_password: str = DB_PASSWORD,
    db_host: str = "localhost",
    db_name: str = "summertest",
) -> Engine:
    """
    Function to create a postgres engine

    :param db_user: User for db
    :param db_password: password for db
    :param db_host: host of db
    :param db_name: name of db
    :return: sqlalchemy engine
    """

    return create_engine(
        f"postgresql+psycopg://{db_user}:{db_password}" f"@/{db_name}",
        future=True,
    )


def create_q3c_extension(__tablename__, ra_column_name, dec_column_name, conn=None):
    print("Executing DDL")

    trig_ddl = DDL(
        "CREATE EXTENSION IF NOT EXISTS q3c;"
        f"CREATE INDEX ON {__tablename__} "
        f"(q3c_ang2ipix({ra_column_name}, {dec_column_name}));"
        f"CLUSTER {__tablename__} USING {__tablename__}_q3c_ang2ipix_idx;"
        f"ANALYZE {__tablename__};"
    )

    if conn is None:
        conn = get_engine().connect()
    conn.execute(trig_ddl)