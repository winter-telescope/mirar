"""
Models for the 'nights' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import REAL, Column, Integer, DATE, Select
from sqlalchemy.orm import relationship, Mapped
from datetime import date

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB, date_field, _exists


class NightsTable(Base):  # pylint: disable=too-few-public-methods
    """
    Nights table in database
    """

    __tablename__ = "nights"
    __table_args__ = {'extend_existing': True}

    nid = Column(Integer, primary_key=True)
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
        return _exists(
            Select(self.sql_model).where(self.sql_model.nightid == self.nightid)
        )
    #
    # def increment_raw(self):
    #     self.rawcount += 1
    #     self.update_entry()


