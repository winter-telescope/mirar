"""
Base class for models
"""

import os

from sqlalchemy.orm import DeclarativeBase

from mirar.database.base_table import BaseTable

DB_NAME = os.getenv("DB_NAME", "winter")


class WinterBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = DB_NAME
