"""
This file contains the credential keys for the database.
"""
import os

DB_USER_KEY = "DB_USER"
DB_PASSWORD_KEY = "DB_PWD"

PG_ADMIN_USER_KEY = "PG_ADMIN_USER"
PG_ADMIN_PWD_KEY = "PG_ADMIN_PWD"

DB_USER = os.getenv(DB_USER_KEY)
DB_PASSWORD = os.getenv(DB_PASSWORD_KEY)

ADMIN_USER = os.getenv(PG_ADMIN_USER_KEY, DB_USER)
ADMIN_PASSWORD = os.getenv(PG_ADMIN_PWD_KEY, DB_PASSWORD)
