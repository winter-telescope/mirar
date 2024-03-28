"""
Module containing postgres util functions
"""

import logging

from sqlalchemy import text
from sqlalchemy_utils import create_database, database_exists

from mirar.database.credentials import (
    DB_HOSTNAME,
    DB_NAME,
    DB_PASSWORD,
    DB_PASSWORD_KEY,
    DB_PORT,
    DB_USER,
    DB_USER_KEY,
)
from mirar.database.engine import get_engine
from mirar.database.errors import DataBaseError

logger = logging.getLogger(__name__)


class PostgresUser:
    """
    Basic Postgres user class for executing functions
    """

    user_env_variable = DB_USER_KEY
    pass_env_variable = DB_PASSWORD_KEY

    def __init__(
        self,
        db_user: str = DB_USER,
        db_password: str = DB_PASSWORD,
        db_hostname: str = DB_HOSTNAME,
        db_name: str = DB_NAME,
        db_port: int = DB_PORT,
    ):
        self.db_user = db_user
        self.db_password = db_password
        self.db_hostname = db_hostname
        self.db_name = db_name
        self.db_port = db_port

    def validate_credentials(self):
        """
        Checks that user credentials exist

        :return: None
        """
        if self.db_user is None:
            err = (
                f"'db_user' is set as None. Please pass a db_user as an argument, "
                f"or set the environment variable '{self.user_env_variable}'."
            )
            logger.error(err)
            raise DataBaseError(err)

        if self.db_password is None:
            err = (
                f"'db_password' is set as None. Please pass a password as an argument, "
                f"or set the environment variable '{self.pass_env_variable}'."
            )
            logger.error(err)
            raise DataBaseError(err)

        engine = get_engine(
            db_name=self.db_name,
            db_user=self.db_user,
            db_password=self.db_password,
            db_hostname=self.db_hostname,
            db_port=self.db_port,
        )

        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))

    @staticmethod
    def create_db(db_name: str):
        """
        Creates a database using credentials, if it does not exist

        :param db_name: DB to create
        :return: None
        """

        engine = get_engine(db_name=db_name)
        if not database_exists(engine.url):
            create_database(engine.url)

        assert database_exists(engine.url)

    @staticmethod
    def has_extension(
        extension_name: str,
        db_name: str,
    ) -> bool:
        """
        Function to check if q3c extension exists

        :param extension_name: name of extension to check
        :param db_name: name of database to check
        :return: boolean extension exists
        """

        engine = get_engine(db_name=db_name)
        with engine.connect() as conn:
            command = text(
                f"SELECT extname FROM pg_extension WHERE extname='{extension_name}'"
            )
            res = conn.execute(command).all()

        assert len(res) <= 1, "More than one extension found"

        return len(res) == 1

    @staticmethod
    def has_schema(
        schema_name: str,
        db_name: str,
    ) -> bool:
        """
        Function to check if schema exists and user can access.

        :param extension_name: name of schema to check
        :param db_name: name of database to check
        :return: boolean extension exists
        """

        engine = get_engine(db_name=db_name)
        with engine.connect() as conn:
            command = text(
                f"SELECT schema_name FROM information_schema.schemata WHERE schema_name='{schema_name}'"
            )
            res = conn.execute(command).all()

        assert len(res) <= 1, "More than one schema found"

        return len(res) == 1
