"""
Base class for models
"""
from sqlalchemy.orm import DeclarativeBase

from mirar.database.base_table import BaseTable

DB_NAME = "winter"


class WinterBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = DB_NAME
