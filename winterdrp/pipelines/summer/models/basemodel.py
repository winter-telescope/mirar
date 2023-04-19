"""
Base class for models
"""
from datetime import date

from pydantic import Field
from sqlalchemy.orm import DeclarativeBase

from winterdrp.processors.sqldatabase.basemodel import BaseTable

db_name = "summer"


class SummerBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = db_name
