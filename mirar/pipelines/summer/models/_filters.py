"""
Models for the 'filters' table
"""
from typing import ClassVar

from pydantic import Field, validator
from sqlalchemy import VARCHAR, Column, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.pipelines.summer.models.base_model import SummerBase
from mirar.processors.sqldatabase.base_model import BaseDB

summer_filters_map = {"u": 1, "g": 2, "r": 3, "i": 4}
fid_field: int = Field(ge=0)


class FiltersTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "filters"

    fuid = Column(Integer, primary_key=True)
    fid = Column(Integer, unique=True)
    filtername = Column(VARCHAR(20), unique=True)
    exposures: Mapped["ExposuresTable"] = relationship(back_populates="filt")


class Filter(BaseDB):
    """
    A pydantic model for a filters database entry
    """

    sql_model: ClassVar = FiltersTable
    fid: int = fid_field
    filtername: str = Field(min_length=1)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(values=self.fid, keys="fid")

    @validator("fid")
    @classmethod
    def validate_fid(cls, field_value):
        """
        Ensure that path exists

        :param field_value: field value
        :return: field value
        """
        assert field_value in list(summer_filters_map.values())
        return field_value


def populate_filters(filter_map: dict = None):
    """
    Iteratively fill up the filters table

    :return: None
    """
    if filter_map is None:
        filter_map = dict(summer_filters_map)

    for filter_name, fid in filter_map.items():
        summer_filter = Filter(fid=fid, filtername=filter_name)
        if not summer_filter.exists():
            summer_filter.insert_entry()
