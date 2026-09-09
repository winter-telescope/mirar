"""
Util functions for database interactions
"""

from sqlalchemy import URL, Engine, NullPool, create_engine

from mirar.database.credentials import DBConfig


def get_engine(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    db_name: str,
    db_user: str | None = None,
    db_password: str | None = None,
    db_hostname: str | None = None,
    db_port: int | None = None,
    db_schema: str | None = None,
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
    config = DBConfig.from_env(
        db_user=db_user,
        db_password=db_password,
        db_hostname=db_hostname,
        db_name=db_name,
        db_port=db_port,
        db_schema=db_schema,
    )

    url_object = URL.create(
        "postgresql+psycopg",
        username=config.db_user,
        password=config.db_password,
        host=config.db_hostname,
        port=config.db_port,
        database=config.db_name,
    )

    return create_engine(
        url_object,
        future=True,
        poolclass=NullPool,
        connect_args={"options": f"-csearch_path={config.db_schema}"},
    )
