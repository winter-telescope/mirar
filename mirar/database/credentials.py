"""
This file contains the credential keys and configuration model for the
database.
"""

import os

from pydantic import BaseModel

DB_USER_KEY = "DB_USER"
DB_PASSWORD_KEY = "DB_PWD"
DB_HOSTNAME_KEY = "DB_HOSTNAME"
DB_NAME_KEY = "DB_NAME"
DB_PORT_KEY = "DB_PORT"
DB_SCHEMA_KEY = "DB_SCHEMA"

PG_ADMIN_USER_KEY = "PG_ADMIN_USER"
PG_ADMIN_PWD_KEY = "PG_ADMIN_PWD"


class DBConfig(BaseModel):
    """
    Database connection settings and credentials.

    Instances should be built with :meth:`from_env`, so that the environment
    is read at call time rather than at import time. This allows credentials
    to be set (e.g. via `os.environ` or a `.env` file) after `mirar` has
    already been imported, and still take effect. Any argument passed to
    `from_env` takes precedence over the corresponding environment variable.
    """

    db_user: str | None = None
    db_password: str | None = None
    db_hostname: str = "127.0.0.1"
    db_name: str = "postgres"
    db_port: int = 5432
    db_schema: str = "public"
    admin_user: str | None = None
    admin_password: str | None = None

    @classmethod
    def from_env(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        cls,
        db_user: str | None = None,
        db_password: str | None = None,
        db_hostname: str | None = None,
        db_name: str | None = None,
        db_port: int | None = None,
        db_schema: str | None = None,
        admin_user: str | None = None,
        admin_password: str | None = None,
    ) -> "DBConfig":
        """
        Build a DBConfig, preferring any argument that is passed in and
        otherwise falling back to the corresponding environment variable.

        :return: DBConfig populated from the arguments/environment
        """
        db_user = db_user if db_user is not None else os.getenv(DB_USER_KEY)
        db_password = (
            db_password if db_password is not None else os.getenv(DB_PASSWORD_KEY)
        )
        return cls(
            db_user=db_user,
            db_password=db_password,
            db_hostname=(
                db_hostname
                if db_hostname is not None
                else os.getenv(DB_HOSTNAME_KEY, "127.0.0.1")
            ),
            db_name=(
                db_name if db_name is not None else os.getenv(DB_NAME_KEY, "postgres")
            ),
            db_port=(
                db_port if db_port is not None else int(os.getenv(DB_PORT_KEY, "5432"))
            ),
            db_schema=(
                db_schema
                if db_schema is not None
                else os.getenv(DB_SCHEMA_KEY, "public")
            ),
            admin_user=(
                admin_user
                if admin_user is not None
                else os.getenv(PG_ADMIN_USER_KEY, db_user)
            ),
            admin_password=(
                admin_password
                if admin_password is not None
                else os.getenv(PG_ADMIN_PWD_KEY, db_password)
            ),
        )
