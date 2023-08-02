"""
Base class for models
"""
from sqlalchemy.orm import DeclarativeBase

from mirar.database.base_model import BaseTable

DB_NAME = "summer"


class SummerBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = DB_NAME
