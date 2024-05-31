"""
Models for the 'exposures' table
"""

import logging
from datetime import date, datetime
from typing import ClassVar

import pandas as pd
from astropy.coordinates import SkyCoord
from pydantic import Field, computed_field
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
from mirar.pipelines.winter.models._fields import FieldsTable, fieldid_field
from mirar.pipelines.winter.models._filters import FiltersTable, fid_field
from mirar.pipelines.winter.models._img_type import ImgTypesTable
from mirar.pipelines.winter.models._nights import Night, NightsTable
from mirar.pipelines.winter.models._programs import Program, default_program
from mirar.pipelines.winter.models.base_model import WinterBase

logger = logging.getLogger(__name__)


class ExposuresTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "exposures"
    __table_args__ = {"extend_existing": True}

    uexpid = Column(
        Integer,
        Sequence(name="exposures_uexpid_seq", start=1, increment=1),
        unique=True,
        autoincrement=True,
        primary_key=True,
    )
    expid = Column(BigInteger, primary_key=False, unique=True, autoincrement=False)
    pipeversion = Column(VARCHAR(10), nullable=True, default=None)
    lastmodified = Column(DateTime(timezone=True))
    # Deterministic ID of exposure

    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))
    filt: Mapped["FiltersTable"] = relationship(back_populates="exposures")

    nightdate: Mapped[int] = mapped_column(ForeignKey("nights.nightdate"))
    night: Mapped["NightsTable"] = relationship(back_populates="exposures")

    fieldid: Mapped[int] = mapped_column(ForeignKey("fields.fieldid"))
    field: Mapped["FieldsTable"] = relationship(back_populates="exposures")

    itid: Mapped[int] = mapped_column(ForeignKey("imgtypes.itid"))
    img_type: Mapped["ImgTypesTable"] = relationship(back_populates="exposures")

    progname: Mapped[str] = mapped_column(ForeignKey("programs.progname"))
    program_name: Mapped["ProgramsTable"] = relationship(back_populates="exposures")

    targname = Column(VARCHAR(MAX_TARGNAME_LEN), nullable=True)
    rawpath = Column(VARCHAR(255), unique=True)

    readoutm = Column(VARCHAR(20), nullable=True)
    readoutv = Column(VARCHAR(10), nullable=True)

    mircover = Column(VARCHAR(10), nullable=True)

    numdiths = Column(Integer, nullable=True)
    dithnum = Column(Integer, nullable=True)
    dithstep = Column(Float, nullable=True)

    galactic_b = Column(Float, nullable=True)
    galactic_l = Column(Float, nullable=True)

    utctime = Column(DateTime(timezone=True))

    exptime = Column(Float, nullable=False)
    expmjd = Column(Float, nullable=False)
    airmass = Column(Float)
    tempture = Column(Float, default=-999)
    windspd = Column(Float, default=-999)
    dewpoint = Column(Float, default=-999)
    humidity = Column(Float, default=-999)
    pressure = Column(Float, default=-999)

    moonaz = Column(Float, default=-999)
    moonalt = Column(Float, default=-999)
    sunalt = Column(Float, default=-999)

    ra = Column(Float)
    dec = Column(Float)
    altitude = Column(Float)
    azimuth = Column(Float)

    ra_column_name = "ra"
    dec_column_name = "dec"

    raw: Mapped["RawsTable"] = relationship(back_populates="exposure_ids")


default_unknown_field = Field(default=-999)


class Exposure(BaseDB):
    """
    A pydantic model for am exposure database entry
    """

    sql_model: ClassVar = ExposuresTable

    expid: int = Field(ge=0)
    fid: int = fid_field
    nightdate: date = Field()  # FIXME : why different to obsdate?
    fieldid: int = fieldid_field
    itid: int = Field(ge=0)
    progname: str = Field(min_length=1)
    targname: str | None = Field(
        min_length=0, max_length=MAX_TARGNAME_LEN, default=None
    )
    rawpath: str = Field(min_length=1)
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
        night = Night(nightdate=self.nightdate)
        logger.debug(f"Searched for night {self.nightdate}")
        if not night.exists():
            night.insert_entry(duplicate_protocol="ignore")

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
