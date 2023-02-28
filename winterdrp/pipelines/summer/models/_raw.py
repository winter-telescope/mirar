"""
Models for the 'raw' table
"""
from typing import ClassVar, List

from pydantic import Field, validator
from datetime import date, datetime
from sqlalchemy import VARCHAR, Column, Integer, ForeignKey
from sqlalchemy.orm import mapped_column, Mapped, relationship
from pathlib import Path

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB, date_field
from winterdrp.pipelines.summer.models._itid import itid_field
from winterdrp.pipelines.summer.models._fields import fieldid_field
from winterdrp.pipelines.summer.models._nights import Nights, NightsTable


class RawTable(Base):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "raw"
    __table_args__ = {'extend_existing': True}

    fid = Column(Integer, primary_key=True)
    filtername = Column(VARCHAR(20), unique=True)
    # nid: Mapped[int] = mapped_column(ForeignKey("nights.nid"))
    nightid: Mapped[int] = mapped_column(ForeignKey("nights.nightid"))
    night: Mapped["NightsTable"] = relationship(back_populates="raw")


class Raw(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = RawTable
    # expid: int = Field(ge=0)
    # savepath: Path = Field(min_length=1)
    # obsdate: date = date_field
    # timeutc: datetime = Field()
    # obsid: int = Field()
    # itid: int = itid_field
    nightid: date = Field() #FIXME : why different to obsdate?
    # fieldid: int = fieldid_field

    # @validator("savepath")
    # @classmethod
    # def validate_time_allocation(cls, field_value: Path):
    #     """
    #     Ensure that path exists
    #
    #     :param field_value: field value
    #     :return: field value
    #     """
    #     assert field_value.exists()
    #     return field_value

    @validator("nightid")
    @classmethod
    def validate_night(cls, field_value: Path):
        """
        Ensure that path exists

        :param field_value: field value
        :return: field value
        """
        night = Nights(nightid=field_value)
        if not night.exists():
            night.insert_entry()
        return field_value
    #
    # def insert_entry(self):
    #     """
    #     Insert the pydantic-ified data into the corresponding sql database
    #
    #     :return: None
    #     """
    #     self._insert_entry()
    #     night = Nights(nightid=self.nightid)
    #     night.increment_raw()


    # def exists(self) -> bool:
    #     """
    #     Checks if the pydantic-ified data exists the corresponding sql database
    #
    #     :return: bool
    #     """
    #     return _exists(Select(self.sql_model).where(self.sql_model.fid == self.fid))
