"""
Module to update database entries
"""
from sqlalchemy import text, update

from mirar.database.base_model import BaseDB
from mirar.database.constraints import DBQueryConstraints
from mirar.database.engine import get_engine


def update_database_entry(
    update_dict: dict,
    db_constraints: DBQueryConstraints,
    db_model: BaseDB,
):
    """
    Update a database entry

    :return: None
    """

    engine = get_engine(db_name=db_model.sql_model.db_name)

    with engine.connect() as conn:
        conn.execute(
            update(db_model.sql_model)
            .where(text(db_constraints.parse_constraints()))
            .values(**update_dict)
        )
        conn.commit()
