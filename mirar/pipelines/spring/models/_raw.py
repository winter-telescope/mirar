"""
Models for the 'raw' table
"""

import logging
import os
from datetime import date, datetime
from typing import ClassVar

import pandas as pd
from astropy.coordinates import SkyCoord
from pydantic import Field, computed_field, field_validator
from sqlalchemy import (  # event,
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
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions import select_from_table
from mirar.paths import __version__

# from mirar.pipelines.spring.models._fields import FieldsTable, fieldid_field
from mirar.pipelines.spring.models._filter import FiltersTable, fid_field
from mirar.pipelines.spring.models._img_type import ImgTypesTable
from mirar.pipelines.spring.models._programs import Program, default_program
from mirar.pipelines.spring.models.base_model import SPRINGBase

logger = logging.getLogger(__name__)
default_unknown_field = Field(default=-999)


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
    img_type: Mapped["ImgTypesTable"] = relationship(back_populates="raws")

    savepath = Column(VARCHAR(255), unique=True)

    ustackid: Mapped[int] = mapped_column(ForeignKey("stacks.ustackid"), nullable=True)
    stacks: Mapped["StacksTable"] = relationship(back_populates="raw")

    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))
    filt: Mapped["FiltersTable"] = relationship(back_populates="raws")

    nightdate = Column(DateTime(timezone=True))
    progname = Column(VARCHAR(255))
    prog_name = relationship("Program", back_populates="raws")
    targname = Column(VARCHAR(MAX_TARGNAME_LEN), nullable=True)
    readoutm = Column(VARCHAR(255), nullable=True)
    readoutv = Column(VARCHAR(255), nullable=True)
    mircover = Column(VARCHAR(255), nullable=True)

    numdiths = Column(Integer, nullable=True)
    dithnum = Column(Integer, nullable=True)
    dithstep = Column(Float, nullable=True, default=None)

    utctime = Column(DateTime(timezone=True))
    exptime = Column(Float, nullable=True)
    expmjd = Column(Float, nullable=True)

    tempture = Column(Float, nullable=True, default=None)
    windspd = Column(Float, nullable=True, default=None)
    dewpoint = Column(Float, nullable=True, default=None)
    humidity = Column(Float, nullable=True, default=None)
    pressure = Column(Float, nullable=True, default=None)

    moonaz = Column(Float, nullable=True, default=None)
    moonalt = Column(Float, nullable=True, default=None)
    sunalt = Column(Float, nullable=True, default=None)

    ra = Column(Float, nullable=True)
    dec = Column(Float, nullable=True)
    altitude = Column(Float, nullable=True)
    azimuth = Column(Float, nullable=True)

    ra_column_name = "ra"
    dec_column_name = "dec"


class Raw(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = RawsTable

    rawid: int = Field(ge=0)
    itid: int = Field(ge=0)
    savepath: str = Field(min_length=1)
    ustackid: int | None = Field(ge=0, default=None)
    fid: int = fid_field
    nightdate: date = Field()  # FIXME : why different to obsdate?
    itid: int = Field(ge=0)
    progname: str = Field(min_length=1)
    targname: str | None = Field(
        min_length=0, max_length=MAX_TARGNAME_LEN, default=None
    )
    readoutm: str | None = Field(default=None)
    readoutv: str | None = Field(default=None)
    mircover: str | None = Field(default=None)

    numdiths: int | None = Field(default=None)
    dithnum: int | None = Field(default=None)
    dithstep: float | None = Field(default=None, ge=0)

    utctime: datetime = Field()
    exptime: float = Field(ge=0)
    expmjd: float = Field(ge=59000)

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

    @computed_field()
    @property
    def galactic_b(self) -> float:
        """
        Returns the galactic latitude of the exposure

        :return: galactic latitude
        """
        return SkyCoord(ra=self.ra, dec=self.dec, unit="deg").galactic.b.deg

    @computed_field()
    @property
    def galactic_l(self) -> float:
        """
        Returns the galactic longitude of the exposure

        :return: galactic longitude
        """
        return SkyCoord(ra=self.ra, dec=self.dec, unit="deg").galactic.l.deg

    @computed_field
    @property
    def pipeversion(self) -> str:
        """
        Returns the version of the pipeline used to process the exposure

        :return: version of the pipeline
        """
        return __version__

    @computed_field
    @property
    def lastmodified(self) -> datetime:
        """
        Returns the current date and time

        :return: current date and time
        """
        return datetime.now()

    def insert_entry(
        self, duplicate_protocol: str, returning_key_names=None
    ) -> pd.DataFrame:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :param duplicate_protocol: protocol to follow if duplicate entry is found
        :param returning_key_names: names of the keys to return
        :return: None
        """
        prog_match = select_from_table(
            DBQueryConstraints(columns="progname", accepted_values=self.progname),
            sql_table=Program.sql_model,
        )
        if prog_match.empty:
            logger.debug(
                f"Program {self.progname} not found in database. "
                f"Using default program {default_program.progname}"
            )
            self.progname = default_program.progname

        return self._insert_entry(
            duplicate_protocol=duplicate_protocol,
            returning_key_names=returning_key_names,
        )

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(values=self.expid, keys="expid")

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
