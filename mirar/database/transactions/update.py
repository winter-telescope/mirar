"""
Module to update database entries
"""
import pandas as pd
from sqlalchemy import column, text, update

from mirar.database.base_table import BaseTable
from mirar.database.constraints import DBQueryConstraints
from mirar.database.engine import get_engine


def _update_database_entry(
    update_dict: dict,
    db_constraints: DBQueryConstraints,
    sql_table: BaseTable,
    returning_key_names: str | list[str] | None = None,
) -> pd.DataFrame:
    """
    Update a database entry

    :return: None
    """

    engine = get_engine(db_name=sql_table.db_name)

    if isinstance(returning_key_names, str):
        returning_key_names = [returning_key_names]

    with engine.connect() as conn:
        stmt = (
            update(sql_table)
            .where(text(db_constraints.parse_constraints()))
            .values(**update_dict)
        )
        if returning_key_names is not None:
            stmt = stmt.returning(*[column(x) for x in returning_key_names])

        res = conn.execute(stmt)
        conn.commit()

        if returning_key_names is not None:
            return pd.DataFrame(res.fetchall())

    return pd.DataFrame()
