"""
This file contains the credential keys and configuration model for the
database.
"""

import os
from typing import Optional

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
    already been imported, and still take effect.
    """

    db_user: Optional[str] = None
    db_password: Optional[str] = None
    db_hostname: str = "127.0.0.1"
    db_name: str = "postgres"
    db_port: int = 5432
    db_schema: str = "public"
    admin_user: Optional[str] = None
    admin_password: Optional[str] = None

    @classmethod
    def from_env(cls) -> "DBConfig":
        """
        Build a DBConfig by reading credentials from the current environment.

        :return: DBConfig populated from environment variables
        """
        db_user = os.getenv(DB_USER_KEY)
        db_password = os.getenv(DB_PASSWORD_KEY)
        return cls(
            db_user=db_user,
            db_password=db_password,
            db_hostname=os.getenv(DB_HOSTNAME_KEY, "127.0.0.1"),
            db_name=os.getenv(DB_NAME_KEY, "postgres"),
            db_port=int(os.getenv(DB_PORT_KEY, 5432)),
            db_schema=os.getenv(DB_SCHEMA_KEY, "public"),
            admin_user=os.getenv(PG_ADMIN_USER_KEY, db_user),
            admin_password=os.getenv(PG_ADMIN_PWD_KEY, db_password),
        )
