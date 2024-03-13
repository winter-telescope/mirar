"""
Module for postgres utilities
"""

import logging

import numpy as np
from sqlalchemy import text

from mirar.database.engine import get_engine

logger = logging.getLogger(__name__)


def get_sequence_key_names_from_table(
    db_table: str,
    db_name: str,
) -> list:
    """
    Gets sequence keys of db table

    :param db_table: database table to use
    :param db_name: database name
    :return: numpy array of keys
    """
    engine = get_engine(db_name=db_name)
    with engine.connect() as conn:
        res = conn.execute(
            text("SELECT c.relname FROM pg_class c WHERE c.relkind = 'S';")
        ).fetchall()

    sequences = [x[0] for x in res]
    seq_tables = np.array([x.split("_")[0] for x in sequences])
    seq_columns = np.array([x.split("_")[1] for x in sequences])
    table_sequence_keys = seq_columns[(seq_tables == db_table)]

    return list(table_sequence_keys)
