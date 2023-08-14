"""
Module to make reference stacks table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, BigInteger, Column, Float, Integer

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.pipelines.winter.models.base_model import WinterBase


class RefStacksTable(WinterBase):
    """
    Table to store Reference Stacks
    """

    __tablename__ = "refstacks"

    stackid = Column(BigInteger, primary_key=True, autoincrement=False)
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
    subdetid = Column(Integer, nullable=True, default=None)
    filter = Column(VARCHAR(10))
    savepath = Column(VARCHAR(255), unique=True)
    coadds = Column(Integer)
    zp = Column(Float)
    zpstd = Column(Float)


class RefStack(BaseDB):
    """
    Pydantic model for a reference stack entry
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
    coadds: int = Field(ge=1)
    zp: float = Field(ge=0)
    zpstd: float = Field(ge=0)

    fieldid: int = Field(ge=0, default=None, nullable=True)
    subdetid: int = Field(ge=0, default=None, nullable=True)
    savepath: str = Field(min_length=1)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(
            values=[self.ra_cent, self.dec_cent], keys=["ra_cent", "dec_cent"]
        )
