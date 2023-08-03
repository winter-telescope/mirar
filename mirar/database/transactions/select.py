"""
Module to select database entries
"""
import pandas as pd
from sqlalchemy import select, text

from mirar.database.base_model import BaseDB
from mirar.database.constraints import DBQueryConstraints
from mirar.database.engine import get_engine


def select_from_table(
    db_constraints: DBQueryConstraints,
    db_model: BaseDB,
    output_columns: list[str] | None = None,
) -> pd.DataFrame:
    """
    Select database entries

    :param db_constraints: database query constraints
    :param db_model: database model
    :param output_columns: columns to output (default: all)
    :return: results
    """

    engine = get_engine(db_name=db_model.sql_model.db_name)

    with engine.connect() as conn:
        stmt = select(db_model.sql_model).where(
            text(db_constraints.parse_constraints())
        )
        res = pd.read_sql(stmt, conn)

    if output_columns is not None:
        res = res[output_columns]

    return res
