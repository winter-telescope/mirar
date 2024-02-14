"""
Models for the 'subdets' table
"""

# pylint: disable=duplicate-code
from typing import ClassVar

from pydantic import Field
from sqlalchemy import Column, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.database.base_model import BaseDB
from mirar.database.transactions import is_populated
from mirar.pipelines.summer.models.base_model import SummerBase

DEFAULT_FIELD = 999999999


class SubdetsTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "subdets"

    qid = Column(Integer, primary_key=True)  # serial counter
    detectorID = Column(Integer)  # WINTER detector ID
    # subdetid = Column(Integer)  # Sub-detector id (1-2)?
    nx = Column(Integer)
    nxtot = Column(Integer)
    ny = Column(Integer)
    nytot = Column(Integer)

    raw: Mapped["RawTable"] = relationship(back_populates="subdets")


class SubDet(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = SubdetsTable
    detectorID: int = Field(Integer, ge=0)
    # subdetid: int = Field(Integer, ge=0)
    nx: int = Field(Integer, ge=0)
    nxtot: int = Field(Integer, ge=0)
    ny: int = Field(Integer, ge=0)
    nytot: int = Field(Integer, ge=0)


def populate_subdets(ndetectors: int = 1, nxtot: int = 1, nytot: int = 1):
    """
    Creates entries in the database based on number of detectors and splits

    :param ndetectors: Number of detectors
    :param nxtot: Number of splits in x direction
    :param nytot: Number of splits in y direction
    :return: None
    """
    if not is_populated(SubdetsTable):
        for ndetector in range(ndetectors):
            for nx in range(nxtot):
                for ny in range(nytot):
                    new = SubDet(
                        detectorID=ndetector,
                        nx=nx + 1,
                        ny=ny + 1,
                        nxtot=nxtot,
                        nytot=nytot,
                    )
                    new.insert_entry()
