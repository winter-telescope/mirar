"""
Models for the 'filters' table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Integer, Select

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB, _exists

summer_filters_map = {"u": 1, "g": 2, "r": 3, "i": 4}


class FiltersTable(Base):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "filters"

    fid = Column(Integer, primary_key=True)
    filtername = Column(VARCHAR(20), unique=True)


class Filters(BaseDB):
    """
    A pydantic model for a filters database entry
    """

    sql_model: ClassVar = FiltersTable
    fid: int = Field(ge=0)
    filtername: str = Field(min_length=1)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return _exists(Select(self.sql_model).where(self.sql_model.fid == self.fid))


def populate_filters(filter_map: dict = None):
    """
    Iteratively fill up the filters table

    :return: None
    """
    if filter_map is None:
        filter_map = dict(summer_filters_map)

    for filter_name, fid in filter_map.items():
        summer_filter = Filters(fid=fid, filtername=filter_name)
        if not summer_filter.exists():
            summer_filter.insert_entry()
