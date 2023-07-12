"""
Module to make reference stacks table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Double, Float, Integer

from mirar.pipelines.winter.models.base_model import WinterBase, dec_field, ra_field
from mirar.processors.sqldatabase.base_model import BaseDB


class RefStacksTable(WinterBase):
    """
    Table to store Reference Stacks
    """

    __tablename__ = "refstacks"

    stackid = Column(Double, primary_key=True, autoincrement=False)
    ra_cent = Column(Float)
    dec_cent = Column(Float)
    ra0_0 = Column(Float)
    dec0_0 = Column(Float)
    ra0_1 = Column(Float)
    dec0_1 = Column(Float)
    ra1_0 = Column(Float)
    dec1_0 = Column(Float)
    ra1_1 = Column(Float)
    dec1_1 = Column(Float)
    fieldid = Column(Integer, nullable=True, default=None)
    subid = Column(Integer, nullable=True, default=None)
    filter = Column(VARCHAR(10))
    savepath = Column(VARCHAR(255))


class RefStacks(BaseDB):
    """
    Pydantic model for reference stacks table
    """

    sql_model: ClassVar = RefStacksTable

    stackid: int = Field(ge=0)
    ra_cent: float = ra_field
    dec_cent: float = dec_field
    ra0_0: float = ra_field
    dec0_0: float = dec_field
    ra0_1: float = ra_field
    dec0_1: float = dec_field
    ra1_0: float = ra_field
    dec1_0: float = dec_field
    ra1_1: float = ra_field
    dec1_1: float = dec_field
    filter: str = Field(min_length=1)

    fieldid: int = Field(ge=0, default=None, nullable=True)
    subid: int = Field(ge=0, default=None, nullable=True)
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
        return self._insert_entry(returning_key_names=returning_keys)
