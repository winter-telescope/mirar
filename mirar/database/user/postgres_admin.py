"""
Postgres Admin class
"""
from sqlalchemy.sql.ddl import DDL

from mirar.database.credentials import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    PG_ADMIN_PWD_KEY,
    PG_ADMIN_USER_KEY,
)
from mirar.database.engine import get_engine
from mirar.database.user.postgres_user import PostgresUser


class PostgresAdmin(PostgresUser):
    """
    An Admin postgres user, with additional functionality for creatying new users
    """

    user_env_varaiable = PG_ADMIN_USER_KEY
    pass_env_variable = PG_ADMIN_PWD_KEY

    def __init__(self, db_user: str = ADMIN_USER, db_password: str = ADMIN_PASSWORD):
        super().__init__(db_user=db_user, db_password=db_password)

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
        )
        with engine.connect() as conn:
            command = DDL(f"CREATE EXTENSION IF NOT EXISTS {extension_name};")
            conn.execute(command)
            conn.commit()

        assert self.has_extension(extension_name=extension_name, db_name=db_name)
