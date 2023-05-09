"""
Models for database and pydantic dataclass models
"""
from typing import Type, Union

from sqlalchemy.orm import DeclarativeBase

from winterdrp.pipelines.reference_building.db_models._ref_components import (
    RefComponents,
    RefComponentsTable,
)
from winterdrp.pipelines.reference_building.db_models._ref_stacks import (
    RefStacks,
    RefStacksTable,
)
from winterdrp.pipelines.reference_building.db_models.basemodel import RefBase
from winterdrp.processors.database.postgres import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_USER,
    PostgresAdmin,
)
from winterdrp.processors.sqldatabase.base_model import BaseTable
from winterdrp.utils.sql import get_engine


def setup_database(base: Union[DeclarativeBase, BaseTable]):
    """
    Function to setup database
    Args:
        base:

    Returns:

    """
    if DB_USER is not None:
        db_name = base.db_name
        admin_engine = get_engine(
            db_name=db_name, db_user=ADMIN_USER, db_password=ADMIN_PASSWORD
        )

        pg_admin = PostgresAdmin()

        if not pg_admin.check_if_db_exists(db_name=db_name):
            pg_admin.create_db(db_name=db_name)

        base.metadata.create_all(
            admin_engine
        )  # extensions need to be created as a superuser

        if not pg_admin.check_if_user_exists(user_name=DB_USER):
            pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)

        pg_admin.grant_privileges(db_name=db_name, db_user=DB_USER)


setup_database(RefBase)
