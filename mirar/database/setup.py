"""
Central function to setup database
"""
import logging
from typing import Union

from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import DeclarativeBase

from mirar.database.base_table import BaseTable
from mirar.database.credentials import DB_PASSWORD, DB_USER
from mirar.database.engine import get_engine
from mirar.database.user import PostgresAdmin, PostgresUser

logger = logging.getLogger(__name__)


def setup_database(db_base: Union[DeclarativeBase, BaseTable]):
    """
    Function to setup database

    :param db_base: BaseTable
    :return: None
    """
    db_name = db_base.db_name
    engine = get_engine(db_name=db_name)

    pg_user = PostgresUser()

    try:
        pg_user.validate_credentials()
    except OperationalError:
        logger.warning(
            f"Failed to credentials for user {DB_USER}. "
            f"Will try creating new user with this name using admin credentials."
        )
        pg_admin = PostgresAdmin()
        pg_admin.validate_credentials()
        pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)
        pg_user.validate_credentials()

    pg_user.create_db(db_name=db_name)

    db_base.metadata.create_all(engine)

    if not pg_user.has_extension(extension_name="q3c", db_name=db_name):
        logger.info(
            f"No 'q3c' extension found. Creating it now for db "
            f"{db_name} with admin user."
        )

        pg_admin = PostgresAdmin()
        pg_admin.validate_credentials()
        pg_admin.create_extension(extension_name="q3c", db_name=db_name)
