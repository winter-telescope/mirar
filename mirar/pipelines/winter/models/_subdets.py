"""
Models for the 'subdets' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import Column, Integer, Select
from sqlalchemy.orm import Mapped, relationship

from mirar.pipelines.winter.constants import subdets
from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.sqldatabase.base_model import BaseDB, _exists
from mirar.utils.sql import get_engine

DEFAULT_FIELD = 999999999


class SubdetsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "subdets"

    subdetid = Column(Integer, primary_key=True)  # serial counter
    boardid = Column(Integer)  # WINTER detector ID
    # subdetid = Column(Integer)  # Sub-detector id (1-2)?
    nx = Column(Integer)
    nxtot = Column(Integer)
    ny = Column(Integer)
    nytot = Column(Integer)

    raw: Mapped["RawTable"] = relationship(back_populates="subdets")


class Subdets(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = SubdetsTable
    boardid: int = Field(Integer, ge=0)
    # subdetid: int = Field(Integer, ge=0)
    nx: int = Field(Integer, ge=1)
    nxtot: int = Field(Integer, ge=1)
    ny: int = Field(Integer, ge=1)
    nytot: int = Field(Integer, ge=1)


def populate_subdets():
    """
    Creates entries in the database based on number of detectors and splits

    """
    engine = get_engine(db_name=SubdetsTable.db_name)
    if not _exists(Select(SubdetsTable), engine=engine):
        for _, row in subdets.iterrows():
            new = Subdets(
                boardid=row["boardid"],
                nx=row["nx"],
                nxtot=row["nxtot"],
                ny=row["ny"],
                nytot=row["nytot"],
            )
            new.insert_entry()
