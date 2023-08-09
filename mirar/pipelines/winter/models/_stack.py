"""
Models for the 'proc' table
"""
import os
from typing import ClassVar

from pydantic import Field, field_validator
from sqlalchemy import REAL, VARCHAR, BigInteger, Column, Integer, Sequence  # event,
from sqlalchemy.orm import relationship

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.pipelines.winter.models.base_model import WinterBase


class StacksTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "stacks"
    __table_args__ = {"extend_existing": True}

    ustackid = Column(
        Integer,
        Sequence(start=1, name="stacks_ustackid_seq"),
        autoincrement=True,
        unique=True,
        primary_key=True,
    )
    stackid = Column(BigInteger, primary_key=False, unique=True, autoincrement=False)
    candidates = relationship("CandidatesTable", back_populates="stack_id")
    raw = relationship("RawsTable", back_populates="stacks")
    diff = relationship("DiffsTable", back_populates="stack_id")
    # procid = Column(Double, primary_key=True, autoincrement=False)

    # rawid: Mapped[int] = Column(Integer, ForeignKey("raw.rawid"), primary_key=True)
    # raw_ids: Mapped["RawTable"] = relationship(back_populates="proc")

    savepath = Column(VARCHAR(255), unique=True)
    wghtpath = Column(VARCHAR(255), unique=True)

    cd1_1 = Column(REAL)
    cd1_2 = Column(REAL)
    cd2_1 = Column(REAL)
    cd2_2 = Column(REAL)
    crval1 = Column(REAL)
    crval2 = Column(REAL)
    crpix1 = Column(REAL)
    crpix2 = Column(REAL)
    zp_auto = Column(REAL)
    fwhm_med = Column(REAL)
    fwhm_std = Column(REAL)
    zp_auto_nstars = Column(Integer)
    zp_auto_std = Column(REAL)
    maglim = Column(REAL)


class Stack(BaseDB):
    """
    A pydantic model for a stack database entry
    """

    sql_model: ClassVar = StacksTable

    # rawid: int = Field(ge=0)
    stackid: int = Field(ge=0)
    savepath: str = Field(min_length=1)
    wghtpath: str = Field(min_length=1)

    cd1_1: float = Field()
    cd1_2: float = Field()
    cd2_1: float = Field()
    cd2_2: float = Field()
    crval1: float = ra_field
    crval2: float = dec_field
    crpix1: float = Field()
    crpix2: float = Field()
    zp_auto: float = Field(ge=0)
    fwhm_med: float = Field(ge=0)
    fwhm_std: float = Field(ge=0)
    zp_auto_nstars: int = Field(ge=0)
    zp_auto_std: float = Field(ge=0)
    maglim: float = Field()

    @field_validator("savepath")
    @classmethod
    def validate_savepath(cls, savepath: str) -> str:
        """
        Ensure that path exists

        :param savepath: savepath
        :return: savepath
        """
        assert os.path.exists(savepath)
        return savepath
