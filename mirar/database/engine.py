"""
Util functions for database interactions
"""

from typing import Optional

from sqlalchemy import URL, Engine, NullPool, create_engine

from mirar.database.credentials import DBConfig


def get_engine(
    db_name: str,
    db_user: Optional[str] = None,
    db_password: Optional[str] = None,
    db_hostname: Optional[str] = None,
    db_port: Optional[int] = None,
    db_schema: Optional[str] = None,
) -> Engine:
    """
    Function to create a postgres engine. Any argument left as None is
    populated from the environment at call time via :class:`DBConfig`.

    :param db_user: User for db
    :param db_password: password for db
    :param db_name: name of db
    :param db_hostname: hostname of db
    :param db_port: port of db
    :param db_schema: schema of db
    :return: sqlalchemy engine
    """
    config = DBConfig.from_env()

    url_object = URL.create(
        "postgresql+psycopg",
        username=db_user if db_user is not None else config.db_user,
        password=db_password if db_password is not None else config.db_password,
        host=db_hostname if db_hostname is not None else config.db_hostname,
        port=db_port if db_port is not None else config.db_port,
        database=db_name,
    )

    schema = db_schema if db_schema is not None else config.db_schema

    return create_engine(
        url_object,
        future=True,
        poolclass=NullPool,
        connect_args={"options": f"-csearch_path={schema}"},
    )
