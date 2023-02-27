"""
Models for the 'field' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import REAL, Column, Integer

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB, dec, ra


class FieldsTable(Base):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "fields"

    fieldid = Column(Integer, primary_key=True)
    fldid = Column(Integer)
    ra = Column(REAL)
    dec = Column(REAL)


class Fields(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = FieldsTable
    fldid: int = Field()
    ra: float = ra
    dec: float = dec
