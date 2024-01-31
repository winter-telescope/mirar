"""
Models for the 'exposures' table
"""

# pylint: disable=duplicate-code
import logging
from datetime import date, datetime
from typing import ClassVar

import pandas as pd
from pydantic import Field, field_validator
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

from mirar.database.base_model import (
    BaseDB,
    alt_field,
    az_field,
    date_field,
    dec_field,
    ra_field,
)
from mirar.pipelines.summer.models._fields import FieldsTable, fieldid_field
from mirar.pipelines.summer.models._filters import FiltersTable, fid_field
from mirar.pipelines.summer.models._img_type import ImgTypesTable
from mirar.pipelines.summer.models._nights import (
    SUMMER_NIGHT_FORMAT,
    Night,
    NightsTable,
)
from mirar.pipelines.summer.models._programs import Program
from mirar.pipelines.summer.models.base_model import SummerBase

logger = logging.getLogger(__name__)


class ExposuresTable(SummerBase):  # pylint: disable=too-few-public-methods
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
    expid = Column(BigInteger, unique=True, autoincrement=False)
    # Deterministic ID of exposure

    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))
    filt: Mapped["FiltersTable"] = relationship(back_populates="exposures")

    nightdate: Mapped[int] = mapped_column(ForeignKey("nights.nightdate"))
    night: Mapped["NightsTable"] = relationship(back_populates="exposures")

    fieldid: Mapped[int] = mapped_column(ForeignKey("fields.fieldid"))
    field: Mapped["FieldsTable"] = relationship(back_populates="exposures")

    itid: Mapped[int] = mapped_column(ForeignKey("imgTypes.itid"))
    img_type: Mapped["ImgTypesTable"] = relationship(back_populates="exposures")

    puid: Mapped[int] = mapped_column(ForeignKey("programs.puid"))
    program_uid: Mapped["ProgramsTable"] = relationship(back_populates="exposures")

    timeutc = Column(DateTime(timezone=True))

    aexptime = Column(Float, nullable=False)
    expmjd = Column(Float, nullable=False)
    airmass = Column(Float)
    shutopen = Column(DateTime(timezone=True))
    shutclsd = Column(DateTime(timezone=True))
    tempture = Column(Float, default=-999)
    windspd = Column(Float, default=-999)
    dewpoint = Column(Float, default=-999)
    humidity = Column(Float, default=-999)
    pressure = Column(Float, default=-999)
    moonra = Column(Float, default=-999)
    moondec = Column(Float, default=-999)
    moonillf = Column(Float, default=-999)
    moonphas = Column(Float, default=-999)
    moonaz = Column(Float, default=-999)
    moonalt = Column(Float, default=-999)
    sunaz = Column(Float, default=-999)
    sunalt = Column(Float, default=-999)
    detsoft = Column(VARCHAR(50), default="unknown")
    detfirm = Column(VARCHAR(50), default="unknown")
    ra = Column(Float)
    dec = Column(Float)
    altitude = Column(Float)
    azimuth = Column(Float)

    ra_column_name = "ra"
    dec_column_name = "dec"

    raw: Mapped["RawTable"] = relationship(back_populates="exposure_ids")
    diff: Mapped["DiffTable"] = relationship(back_populates="exposure_ids")


default_unknown_field = Field(default=-999)


class Exposure(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = ExposuresTable

    expid: int = Field(ge=0)
    fid: int = fid_field
    nightdate: date = date_field
    fieldid: int = fieldid_field
    itid: int = Field(ge=0)
    puid: int = Field(ge=0)
    timeutc: datetime = Field()
    aexptime: float = Field(ge=0)
    expmjd: float = Field(ge=59000)
    airmass: float = Field(ge=1.0)
    shutopen: datetime = Field()
    shutclsd: datetime = Field()
    tempture: float = default_unknown_field
    windspd: float = default_unknown_field
    dewpoint: float = default_unknown_field
    humidity: float = default_unknown_field
    pressure: float = default_unknown_field
    moonra: float = Field(ge=0.0, le=360.0, default=None)
    moondec: float = Field(title="Dec (degrees)", ge=-90.0, le=90, default=None)
    moonillf: float = default_unknown_field
    moonphas: float = default_unknown_field
    moonaz: float = default_unknown_field
    moonalt: float = default_unknown_field
    sunaz: float = default_unknown_field
    sunalt: float = default_unknown_field
    detfirm: str = Field(default="unknown")
    detsoft: str = Field(default="unknown")
    ra: float = ra_field
    dec: float = dec_field
    altitude: float = alt_field
    azimuth: float = az_field

    @field_validator("nightdate", mode="before")
    @classmethod
    def validate_fid(cls, nightdate: str) -> datetime:
        """
        Ensure that path exists

        :param nightdate: str nightdate
        :return: datetime nightdate
        """
        return datetime.strptime(nightdate, SUMMER_NIGHT_FORMAT)

    def insert_entry(self, returning_key_names=None) -> pd.DataFrame:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: None
        """
        night = Night(nightdate=self.nightdate)
        if not night.exists():
            night.insert_entry()

        logger.debug(f"puid: {self.puid}")
        if not Program._exists(values=self.puid, keys="puid"):
            self.puid = 1

        return self._insert_entry()

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(values=self.expid, keys="expid")
