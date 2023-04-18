"""
Base class for models
"""
from sqlalchemy.orm import DeclarativeBase

from winterdrp.processors.sqldatabase.basemodel import BaseTable

db_name = "winter"


class WinterBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = db_name
