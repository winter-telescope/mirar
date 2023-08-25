"""
Util functions for database interactions
"""
from sqlalchemy import Engine, NullPool, create_engine

from mirar.database.credentials import DB_HOSTNAME, DB_PASSWORD, DB_PORT, DB_USER


def get_engine(
    db_name: str,
    db_user: str = DB_USER,
    db_password: str = DB_PASSWORD,
    db_hostname: str = DB_HOSTNAME,
    db_port: int = DB_PORT,
) -> Engine:
    """
    Function to create a postgres engine

    :param db_user: User for db
    :param db_password: password for db
    :param db_name: name of db
    :param db_hostname: hostname of db
    :param db_port: port of db
    :return: sqlalchemy engine
    """

    return create_engine(
        f"postgresql+psycopg://{db_user}:{db_password}@{db_hostname}:{db_port}/{db_name}",
        future=True,
        echo=True,
        poolclass=NullPool,
    )
