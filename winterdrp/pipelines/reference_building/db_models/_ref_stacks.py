"""
Module to make reference stacks table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Float, Integer

from winterdrp.pipelines.reference_building.db_models.basemodel import (
    RefBase,
    dec_field,
    ra_field,
)
from winterdrp.processors.sqldatabase.basemodel import BaseDB


class RefStacksTable(RefBase):
    """
    Table to store Reference Stacks
    """

    __tablename__ = "refstacks"

    stackid = Column(Integer, primary_key=True)
    ra_cent = Column(Float)
    dec_cent = Column(Float)
    savepath = Column(VARCHAR(255))


class RefStacks(BaseDB):
    """
    Pydantic model for reference stacks table
    """

    sql_model: ClassVar = RefStacksTable

    ra_cent: float = ra_field
    dec_cent: float = dec_field
    savepath: str = Field(min_length=1)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(
            values=[self.ra_cent, self.dec_cent], keys=["ra_cent", "dec_cent"]
        )

    def insert_entry(self, returning_keys=None):
        return self._insert_entry(returning_keys=returning_keys)
