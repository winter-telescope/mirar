"""
Models for the 'imgTypes' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.sqldatabase.base_model import BaseDB

ALL_ITID = ["SCIENCE", "CAL", "FOCUS", "POINTING", "NULL"]
itid_dict = {}
for i, img_type in enumerate(ALL_ITID):
    itid_dict[img_type] = i + 1
DEFAULT_ITID = 5

itid_field = Field(ge=0, le=5, default=5)


class ImgTypesTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    ITID table in database
    """

    __tablename__ = "imgTypes"

    itid = Column(Integer, primary_key=True)
    imgtype = Column(VARCHAR(20), unique=True)
    exposures: Mapped["ExposuresTable"] = relationship(back_populates="img_type")


class ImgTypes(BaseDB):
    """
    A pydantic model for a ITID database entry
    """

    sql_model: ClassVar = ImgTypesTable
    itid: int = itid_field
    imgtype: str = Field(min_length=1, max_length=20)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(values=self.itid, keys="itid")


def populate_itid():
    """
    Iteratively fill up the ITID table

    :return: None
    """
    for ind, imgtype in enumerate(ALL_ITID):
        itid = ImgTypes(itid=ind + 1, imgtype=imgtype)
        if not itid.exists():
            itid.insert_entry()
