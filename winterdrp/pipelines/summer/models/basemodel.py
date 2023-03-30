"""
Base class for models
"""
from datetime import date
from typing import ClassVar

from pydantic import BaseModel, Extra, Field, root_validator, validator
from sqlalchemy import Insert, Select, Table, inspect, select
from sqlalchemy.orm import DeclarativeBase

from winterdrp.processors.sqldatabase.basemodel import BaseTable
from winterdrp.utils.sql import get_engine

db_name = "summertest"


class SummerBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = db_name


ra_field: float = Field(title="RA (degrees)", ge=0.0, le=360.0)
dec_field: float = Field(title="Dec (degrees)", ge=-90.0, le=90.0)
alt_field: float = Field(title="Alt (degrees)", ge=0.0, le=90.0)
az_field: float = Field(title="Az (degrees)", ge=0.0, le=360.0)
date_field: date = Field()
