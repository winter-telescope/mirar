"""
Models for the 'nights' table
"""
from datetime import date
from typing import ClassVar

from sqlalchemy import DATE, Column, Integer, Sequence
from sqlalchemy.orm import Mapped, relationship

from winterdrp.pipelines.summer.models.basemodel import SummerBase, date_field
from winterdrp.processors.sqldatabase.basemodel import BaseDB


class NightsTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Nights table in database
    """

    __tablename__ = "nights"
    __table_args__ = {"extend_existing": True}

    nid = Column(Integer, primary_key=True)  # Starts at 1 and
    # creates a sequence by default
    nightid = Column(DATE, unique=True)
    raw: Mapped["RawTable"] = relationship(back_populates="night")


class Nights(BaseDB):
    """
    A pydantic model for a nights database entry
    """

    sql_model: ClassVar = NightsTable
    nightid: date = date_field

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(value=self.nightid, key="nightid")

    #
    # def increment_raw(self):
    #     self.rawcount += 1
    #     self.update_entry()
