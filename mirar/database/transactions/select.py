"""
Module to select database entries
"""
from sqlalchemy import select, text
from sqlalchemy.orm import load_only

from mirar.database.base_model import BaseDB
from mirar.database.constraints import DBQueryConstraints
from mirar.database.engine import get_engine


def select_from_table(
    db_constraints: DBQueryConstraints,
    db_model: BaseDB,
    output_columns: list[str] = None,
) -> list[tuple]:
    """
    Select database entries

    :return: results
    """

    engine = get_engine(db_name=db_model.sql_model.db_name)

    with engine.connect() as conn:
        if output_columns is not None:
            stmt = select(db_model.sql_model).options(load_only(*output_columns))
        else:
            stmt = select(db_model.sql_model)

        res = conn.execute(
            stmt.where(text(db_constraints.parse_constraints()))
        ).fetchall()
    return res
