"""
Module for postgres utilities
"""
import logging
import os

import numpy as np
import psycopg

logger = logging.getLogger(__name__)


DB_USER_KEY = "DB_USER"
DB_PASSWORD_KEY = "DB_PWD"

PG_ADMIN_USER_KEY = "PG_ADMIN_USER"
PG_ADMIN_PWD_KEY = "PG_ADMIN_PWD"

DB_USER = os.getenv(DB_USER_KEY)
DB_PASSWORD = os.getenv(DB_PASSWORD_KEY)

ADMIN_USER = os.getenv(PG_ADMIN_USER_KEY, DB_USER)
ADMIN_PASSWORD = os.getenv(PG_ADMIN_PWD_KEY, DB_PASSWORD)

POSTGRES_DUPLICATE_PROTOCOLS = ["fail", "ignore", "replace"]


def get_sequence_key_names_from_table(
    db_table: str,
    db_name: str,
    db_user: str = DB_USER,
    db_password: str = DB_PASSWORD,
) -> np.ndarray:
    """
    Gets sequence keys of db table

    :param db_table: database table to use
    :param db_name: database name
    :param db_user: database user
    :param db_password: dataname password
    :return: numpy array of keys

    """
    with psycopg.connect(
        f"dbname={db_name} user={db_user} password={db_password}"
    ) as conn:
        conn.autocommit = True
        sequences = [
            x[0]
            for x in conn.execute(
                "SELECT c.relname FROM pg_class c WHERE c.relkind = 'S';"
            ).fetchall()
        ]
        seq_tables = np.array([x.split("_")[0] for x in sequences])
        seq_columns = np.array([x.split("_")[1] for x in sequences])
        table_sequence_keys = seq_columns[(seq_tables == db_table)]
    return table_sequence_keys
