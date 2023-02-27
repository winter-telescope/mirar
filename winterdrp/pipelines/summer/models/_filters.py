"""
Models for the 'filters' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Integer

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB


class FiltersTable(Base):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "fields"

    fid = Column(Integer, primary_key=True)
    filterid = Column(Integer)
    filtername = Column(VARCHAR(20))


class Filters(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = FiltersTable
    filterid: int = Field(ge=0)
    filtername: str = Field(min_length=1)
