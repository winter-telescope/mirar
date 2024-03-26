"""
This file contains the credential keys for the database.
"""

import os

DB_USER_KEY = "DB_USER"
DB_PASSWORD_KEY = "DB_PWD"
DB_HOSTNAME_KEY = "DB_HOSTNAME"
DB_NAME_KEY = "DB_NAME"
DB_PORT_KEY = "DB_PORT"
DB_SCHEMA_KEY = "DB_SCHEMA"

PG_ADMIN_USER_KEY = "PG_ADMIN_USER"
PG_ADMIN_PWD_KEY = "PG_ADMIN_PWD"

DB_USER = os.environ.get(DB_USER_KEY)
DB_PASSWORD = os.environ.get(DB_PASSWORD_KEY)
DB_HOSTNAME = os.environ.get(DB_HOSTNAME_KEY, "127.0.0.1")
DB_NAME = os.environ.get(DB_NAME_KEY, "postgres")
DB_PORT = os.environ.get(DB_PORT_KEY, 5432)
DB_SCHEMA = os.environ.get(DB_SCHEMA_KEY, "public")

ADMIN_USER = os.environ.get(PG_ADMIN_USER_KEY, DB_USER)
ADMIN_PASSWORD = os.environ.get(PG_ADMIN_PWD_KEY, DB_PASSWORD)
