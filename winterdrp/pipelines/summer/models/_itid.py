"""
Models for the 'itid' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Integer, Select
from sqlalchemy.orm import Mapped, relationship

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB, _exists

ALL_ITID = ["SCIENCE", "CAL", "FOCUS", "POINTING", "OTHER"]
DEFAULT_ITID = 5

itid_field = Field(ge=0, le=5, default=5)


class ITIDsTable(Base):  # pylint: disable=too-few-public-methods
    """
    ITID table in database
    """

    __tablename__ = "itid"

    itid = Column(Integer, primary_key=True)
    imgtype = Column(VARCHAR(20), unique=True)
    raw: Mapped["RawTable"] = relationship(back_populates="img_type")


class ITIDs(BaseDB):
    """
    A pydantic model for a ITID database entry
    """

    sql_model: ClassVar = ITIDsTable
    itid: int = itid_field
    imgtype: str = Field(min_length=1, max_length=20)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return _exists(Select(self.sql_model).where(self.sql_model.itid == self.itid))


def populate_itid():
    """
    Iteratively fill up the ITID table

    :return: None
    """
    for i, img_type in enumerate(ALL_ITID):
        itid = ITIDs(itid=i + 1, imgtype=img_type)
        if not itid.exists():
            itid.insert_entry()
