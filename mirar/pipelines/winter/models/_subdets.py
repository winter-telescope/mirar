"""
Models for the 'subdets' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import Column, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.database.base_model import BaseDB
from mirar.database.transactions import is_populated
from mirar.pipelines.winter.constants import subdets
from mirar.pipelines.winter.models.base_model import WinterBase


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

    raw: Mapped["RawsTable"] = relationship(back_populates="subdets")


class Subdet(BaseDB):
    """
    A pydantic model for a subdet database entry
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
    if not is_populated(SubdetsTable):
        for _, row in subdets.iterrows():
            new = Subdet(
                boardid=row["boardid"],
                nx=row["nx"],
                nxtot=row["nxtot"],
                ny=row["ny"],
                nytot=row["nytot"],
            )
            new.insert_entry()
