"""
Models for the 'raw' table
"""

import os
from datetime import datetime
from typing import ClassVar

from astropy.coordinates import SkyCoord
from pydantic import Field, computed_field, field_validator
from sqlalchemy import (
    VARCHAR,
    BigInteger,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    Sequence,
)
from sqlalchemy.orm import Mapped, mapped_column, relationship
from wintertoo.data import MAX_TARGNAME_LEN

from mirar.database.base_model import BaseDB, alt_field, az_field, dec_field, ra_field
from mirar.paths import __version__
from mirar.pipelines.spring.models._exposures import default_unknown_field
from mirar.pipelines.spring.models._filters import FiltersTable, fid_field
from mirar.pipelines.spring.models.base_model import SPRINGBase


class RawsTable(SPRINGBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "raws"
    __table_args__ = {"extend_existing": True}

    urawid = Column(
        Integer,
        Sequence(start=1, name="raw_urawid_seq"),
        autoincrement=True,
        primary_key=True,
    )
    rawid = Column(BigInteger, primary_key=False, unique=True, autoincrement=False)

    itid: Mapped[int] = mapped_column(ForeignKey("imgtypes.itid"))
    img_type: Mapped["ImgTypesTable"] = relationship()

    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))
    filt: Mapped["FiltersTable"] = relationship()

    progname: Mapped[str] = mapped_column(ForeignKey("programs.progname"))
    program_name: Mapped["ProgramsTable"] = relationship()

    savepath = Column(VARCHAR(255), unique=True)

    ustackid: Mapped[int] = mapped_column(ForeignKey("stacks.ustackid"), nullable=True)
    stacks: Mapped["StacksTable"] = relationship(back_populates="raws")

    targname = Column(VARCHAR(MAX_TARGNAME_LEN), nullable=True)
    readoutm = Column(VARCHAR(255), nullable=True)
    readoutv = Column(VARCHAR(255), nullable=True)
    mircover = Column(VARCHAR(255), nullable=True)

    numdiths = Column(Integer, nullable=True)
    dithnum = Column(Integer, nullable=True)
    dithstep = Column(Float, nullable=True)

    utctime = Column(DateTime(timezone=True), nullable=True)
    exptime = Column(Float, nullable=True)
    expmjd = Column(Float, nullable=True)

    tempture = Column(Float, default=-999)
    windspd = Column(Float, default=-999)
    dewpoint = Column(Float, default=-999)
    humidity = Column(Float, default=-999)
    pressure = Column(Float, default=-999)

    moonaz = Column(Float, default=-999)
    moonalt = Column(Float, default=-999)
    sunalt = Column(Float, default=-999)

    galactic_b = Column(Float, nullable=True)
    galactic_l = Column(Float, nullable=True)

    pipeversion = Column(VARCHAR(10), nullable=True, default=None)
    lastmodified = Column(DateTime(timezone=True), nullable=True)

    ra = Column(Float)
    dec = Column(Float)
    altitude = Column(Float)
    azimuth = Column(Float)

    ra_column_name = "ra"
    dec_column_name = "dec"


class Raw(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = RawsTable

    rawid: int = Field(ge=0)
    itid: int = Field(ge=0)
    fid: int = fid_field
    progname: str = Field(min_length=1)
    savepath: str = Field(min_length=1)

    targname: str | None = Field(
        min_length=0, max_length=MAX_TARGNAME_LEN, default=None
    )
    ustackid: int | None = Field(ge=0, default=None)

    readoutm: str | None = Field(default=None)
    readoutv: str | None = Field(default=None)
    mircover: str | None = Field(default=None)

    numdiths: int | None = Field(default=None)
    dithnum: int | None = Field(default=None)
    dithstep: float | None = Field(default=None, ge=0)

    utctime: datetime | None = Field(default=None)
    exptime: float | None = Field(default=None, ge=0)
    expmjd: float | None = Field(default=None)

    tempture: float = default_unknown_field
    windspd: float = default_unknown_field
    dewpoint: float = default_unknown_field
    humidity: float = default_unknown_field
    pressure: float = default_unknown_field

    moonaz: float = default_unknown_field
    moonalt: float = default_unknown_field
    sunalt: float = default_unknown_field

    ra: float = ra_field
    dec: float = dec_field
    altitude: float = alt_field
    azimuth: float = az_field

    @computed_field
    @property
    def galactic_b(self) -> float:
        """
        Returns the galactic latitude of the raw image pointing
        """
        return SkyCoord(ra=self.ra, dec=self.dec, unit="deg").galactic.b.deg

    @computed_field
    @property
    def galactic_l(self) -> float:
        """
        Returns the galactic longitude of the raw image pointing
        """
        return SkyCoord(ra=self.ra, dec=self.dec, unit="deg").galactic.l.deg

    @computed_field
    @property
    def pipeversion(self) -> str:
        """
        Returns the version of the pipeline used
        """
        return __version__

    @computed_field
    @property
    def lastmodified(self) -> datetime:
        """
        Returns the current date and time
        """
        return datetime.now()

    @field_validator("savepath")
    @classmethod
    def validate_savepath(cls, savepath: str) -> str:
        """
        Ensure that path exists
        """
        assert os.path.exists(savepath)
        return savepath

    # @field_validator("uexpid")
    # @classmethod
    # def validate_expid(cls, uexpid: int) -> int:
    #     """
    #     Ensure that expid exists in exposures table

    #     :param uexpid: unique exposure id
    #     """
    #     assert Exposure._exists(keys="uexpid", values=uexpid)
    #     return uexpid
