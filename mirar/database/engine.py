"""
Util functions for database interactions
"""

from sqlalchemy import URL, Engine, NullPool, create_engine

from mirar.database.credentials import (
    DB_HOSTNAME,
    DB_PASSWORD,
    DB_PORT,
    DB_SCHEMA,
    DB_USER,
)


def get_engine(
    db_name: str,
    db_user: str = DB_USER,
    db_password: str = DB_PASSWORD,
    db_hostname: str = DB_HOSTNAME,
    db_port: int = DB_PORT,
    db_schema: str = DB_SCHEMA,
) -> Engine:
    """
    Function to create a postgres engine

    :param db_user: User for db
    :param db_password: password for db
    :param db_name: name of db
    :param db_hostname: hostname of db
    :param db_port: port of db
    :param db_schema: schema of db
    :return: sqlalchemy engine
    """

    url_object = URL.create(
        "postgresql+psycopg",
        username=db_user,
        password=db_password,
        host=db_hostname,
        port=db_port,
        database=db_name,
    )

    return create_engine(
        url_object,
        future=True,
        poolclass=NullPool,
        connect_args={"options": f"-csearch_path={db_schema}"},
    )
