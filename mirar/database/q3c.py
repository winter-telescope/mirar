"""
Module to create q3c extension and index on table
"""
import logging

from sqlalchemy import DDL

from mirar.database.engine import get_engine

logger = logging.getLogger(__name__)


def create_q3c_extension(
    db_name: str,
    table_name: str,
    ra_column_name: str,
    dec_column_name: str,
):
    """
    Function to create q3c extension and index on table

    :param db_name: Name of database
    :param table_name: Name of table
    :param ra_column_name: ra column name
    :param dec_column_name: dec column name
    :return:
    """

    logger.info(f"Creating q3c extension and index on table {table_name}...")

    trig_ddl = DDL(
        f"CREATE INDEX ON {table_name} "
        f"(q3c_ang2ipix({ra_column_name}, {dec_column_name}));"
        f"CLUSTER {table_name} USING {table_name}_q3c_ang2ipix_idx;"
        f"ANALYZE {table_name};"
    )

    engine = get_engine(db_name=db_name)

    with engine.connect() as conn:
        conn.execute(trig_ddl)
        conn.commit()
