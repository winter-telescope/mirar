"""
Module to update database entries
"""

from sqlalchemy import text, update

from mirar.database.base_table import BaseTable
from mirar.database.constraints import DBQueryConstraints
from mirar.database.engine import get_engine


def _update_database_entry(
    update_dict: dict,
    db_constraints: DBQueryConstraints,
    sql_table: BaseTable,
):
    """
    Update a database entry

    :param update_dict: dictionary of values to update
    :param db_constraints: database query constraints
    :param sql_table: database SQL table
    :return: None
    """

    engine = get_engine(db_name=sql_table.db_name)

    with engine.connect() as conn:
        stmt = (
            update(sql_table)
            .where(text(db_constraints.parse_constraints()))
            .values(**update_dict)
        )

        conn.execute(stmt)
        conn.commit()
