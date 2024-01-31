"""
Models for the 'filters' table
"""

# pylint: disable=duplicate-code
from typing import ClassVar

from pydantic import Field, field_validator
from sqlalchemy import VARCHAR, Column, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.database.base_model import BaseDB
from mirar.pipelines.summer.models.base_model import SummerBase

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
        return self._exists(values=self.fid, keys="fid")

    @field_validator("fid")
    @classmethod
    def validate_fid(cls, filter_id):
        """
        Ensure that path exists

        :param filter_id: field value
        :return: field value
        """
        assert filter_id in list(summer_filters_map.values())
        return filter_id


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
