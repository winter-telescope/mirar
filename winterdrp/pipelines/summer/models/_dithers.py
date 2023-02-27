"""
Models for the 'dithers' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import REAL, Column, Integer

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB


class DithersTable(Base):  # pylint: disable=too-few-public-methods
    """
    Dithers table in database
    """

    __tablename__ = "dithers"

    did = Column(Integer, primary_key=True)
    deltara = Column(REAL)
    deltadec = Column(REAL)


class Dithers(BaseDB):
    """
    A pydantic model for a dithers database entry
    """

    sql_model: ClassVar = DithersTable
    deltara: float = Field()
    deltadec: float = Field()
