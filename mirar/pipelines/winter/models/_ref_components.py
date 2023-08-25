"""
Module to make reference components table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, BigInteger, Column, Float, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.pipelines.winter.models.base_model import WinterBase


class RefComponentsTable(WinterBase):
    """
    Table for individual reference images
    """

    __tablename__ = "refcomponents"

    compid = Column(BigInteger, primary_key=True, autoincrement=False)
    queries: Mapped["RefQueriesTable"] = relationship(back_populates="components")

    mfid = Column(BigInteger)
    xtnsnid = Column(Integer)
    ra_cent = Column(Float)
    dec_cent = Column(Float)
    ramin = Column(Float)
    ramax = Column(Float)
    decmin = Column(Float)
    decmax = Column(Float)
    savepath = Column(VARCHAR(255), unique=True)
    filter = Column(VARCHAR(10))
    zp = Column(Float)
    zpstd = Column(Float)
    seeing = Column(Float)


class RefComponent(BaseDB):
    """
    Pydantic model for a Reference component entry
    """

    sql_model: ClassVar = RefComponentsTable
    compid: int = Field(ge=0)
    savepath: str = Field(min_length=1)
    mfid: int = Field(ge=0)
    xtnsnid: int = Field(ge=0)
    ramin: float = ra_field
    ramax: float = ra_field
    decmin: float = dec_field
    decmax: float = dec_field
    ra_cent: float = ra_field
    dec_cent: float = dec_field
    filter: str = Field(min_length=1)
    zp: float = Field(ge=0)
    zpstd: float = Field(ge=0)
    seeing: float = Field()

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(values=[self.mfid, self.xtnsnid], keys=["mfid", "xtnsnid"])
