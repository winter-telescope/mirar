"""
Util functions for database interactions
"""
from sqlalchemy import DDL, Engine, NullPool, create_engine

from mirar.processors.database.postgres import DB_PASSWORD, DB_USER


def get_engine(
    db_user: str = DB_USER,
    db_password: str = DB_PASSWORD,
    db_name: str = "summer",
) -> Engine:
    """
    Function to create a postgres engine

    :param db_user: User for db
    :param db_password: password for db
    :param db_name: name of db
    :return: sqlalchemy engine
    """

    return create_engine(
        f"postgresql+psycopg://{db_user}:{db_password}" f"@/{db_name}",
        future=True,
        poolclass=NullPool,
    )


def create_q3c_extension(
    table_name: str, ra_column_name: str, dec_column_name: str, conn=None
):
    """
    Function to create q3c extension and index on table

    :param table_name: Name of table
    :param ra_column_name: ra column name
    :param dec_column_name: dec column name
    :param conn: connection to db
    :return:
    """

    trig_ddl = DDL(
        "CREATE EXTENSION IF NOT EXISTS q3c;"
        f"CREATE INDEX ON {table_name} "
        f"(q3c_ang2ipix({ra_column_name}, {dec_column_name}));"
        f"CLUSTER {table_name} USING {table_name}_q3c_ang2ipix_idx;"
        f"ANALYZE {table_name};"
    )

    if conn is None:
        engine = get_engine()
        with engine.connect() as new_conn:
            new_conn.execute(trig_ddl)
            new_conn.commit()
    else:
        conn.execute(trig_ddl)
        conn.commit()
