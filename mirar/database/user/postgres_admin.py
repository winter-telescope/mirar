"""
Postgres Admin class
"""

from sqlalchemy.sql.ddl import DDL

from mirar.database.credentials import PG_ADMIN_PWD_KEY, PG_ADMIN_USER_KEY, DBConfig
from mirar.database.engine import get_engine
from mirar.database.user.postgres_user import PostgresUser


class PostgresAdmin(PostgresUser):
    """
    An Admin postgres user, with additional functionality for creating new users
    """

    user_env_variable = PG_ADMIN_USER_KEY
    pass_env_variable = PG_ADMIN_PWD_KEY

    def __init__(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        db_user: str | None = None,
        db_password: str | None = None,
        db_hostname: str | None = None,
        db_name: str | None = None,
        db_port: int | None = None,
    ):
        config = DBConfig.from_env(
            admin_user=db_user,
            admin_password=db_password,
            db_hostname=db_hostname,
            db_name=db_name,
            db_port=db_port,
        )
        super().__init__(
            db_user=config.admin_user,
            db_password=config.admin_password,
            db_hostname=config.db_hostname,
            db_name=config.db_name,
            db_port=config.db_port,
        )

    def create_new_user(self, new_db_user: str, new_password: str):
        """
        Create a new postgres user

        :param new_db_user: new username
        :param new_password: new user password
        :return: None
        """
        engine = get_engine(
            db_name="postgres",
            db_user=self.db_user,
            db_password=self.db_password,
            db_hostname=self.db_hostname,
            db_port=self.db_port,
        )
        with engine.connect() as conn:
            command = DDL(
                f"CREATE ROLE {new_db_user} WITH password '{new_password}' CREATEDB NOCREATEROLE LOGIN;"
            )
            conn.execute(command)
            conn.commit()

    def create_extension(self, extension_name: str, db_name: str):
        """
        Function to create new extension for database

        :param extension_name: name of extension to create
        :param db_name: name of database to create extension in
        :return: None
        """
        engine = get_engine(
            db_name=db_name,
            db_user=self.db_user,
            db_password=self.db_password,
            db_hostname=self.db_hostname,
            db_port=self.db_port,
        )
        with engine.connect() as conn:
            command = DDL(f"CREATE EXTENSION IF NOT EXISTS {extension_name};")
            conn.execute(command)
            conn.commit()

        assert self.has_extension(extension_name=extension_name, db_name=db_name)

    def create_schema(self, schema_name: str, db_name: str, db_user: str):
        """
        Function to create new schema for database

        :param extension_name: name of schema to create
        :param db_name: name of database to create schema in
        :param db_user: name of schema owner
        :return: None
        """
        engine = get_engine(
            db_name=db_name,
            db_user=self.db_user,
            db_password=self.db_password,
            db_hostname=self.db_hostname,
            db_port=self.db_port,
        )
        with engine.connect() as conn:
            command = DDL(
                f"CREATE SCHEMA IF NOT EXISTS {schema_name} AUTHORIZATION {db_user};"
            )
            conn.execute(command)
            conn.commit()

        assert self.has_schema(schema_name=schema_name, db_name=db_name)
