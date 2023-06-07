"""
Models for the 'exposures' table
"""
# from mirar.utils.sql import create_q3c_extension
import logging
from datetime import date, datetime
from typing import ClassVar

from pydantic import Field
from sqlalchemy import (  # event,
    VARCHAR,
    Column,
    DateTime,
    Double,
    Float,
    ForeignKey,
    Integer,
    Sequence,
)
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.pipelines.summer.models._fields import FieldsTable, fieldid_field
from mirar.pipelines.summer.models._filters import FiltersTable, fid_field
from mirar.pipelines.summer.models._img_type import ImgTypesTable
from mirar.pipelines.summer.models._nights import Night, NightsTable
from mirar.pipelines.summer.models._programs import ProgramsTable, default_program
from mirar.pipelines.summer.models.base_model import SummerBase
from mirar.processors.sqldatabase.base_model import (
    BaseDB,
    alt_field,
    az_field,
    dec_field,
    ra_field,
)

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
    )
    expid = Column(Double, primary_key=True, unique=True, autoincrement=False)
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

    AExpTime = Column(Float, nullable=False)
    expMJD = Column(Float, nullable=False)
    airmass = Column(Float)
    shutopen = Column(DateTime(timezone=True))
    shutclsd = Column(DateTime(timezone=True))
    tempture = Column(Float, default=-999)
    windspd = Column(Float, default=-999)
    Dewpoint = Column(Float, default=-999)
    Humidity = Column(Float, default=-999)
    Pressure = Column(Float, default=-999)
    Moonra = Column(Float, default=-999)
    Moondec = Column(Float, default=-999)
    Moonillf = Column(Float, default=-999)
    Moonphas = Column(Float, default=-999)
    Moonaz = Column(Float, default=-999)
    Moonalt = Column(Float, default=-999)
    Sunaz = Column(Float, default=-999)
    Sunalt = Column(Float, default=-999)
    Detsoft = Column(VARCHAR(50), default="unknown")
    Detfirm = Column(VARCHAR(50), default="unknown")
    ra = Column(Float)
    dec = Column(Float)
    altitude = Column(Float)
    azimuth = Column(Float)

    ra_column_name = "ra"
    dec_column_name = "dec"

    raw: Mapped["RawTable"] = relationship(back_populates="exposure_ids")


# @event.listens_for(target=RawTable.__table__, identifier="after_create")
# def raw_q3c(tbl, conn, *args, **kw):
#     create_q3c_extension(
#         conn=conn,
#         __tablename__=RawTable.__tablename__,
#         ra_column_name=RawTable.ra_column_name,
#         dec_column_name=RawTable.dec_column_name,
#     )


default_unknown_field = Field(default=-999)


class Exposure(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = ExposuresTable

    expid: int = Field(ge=0)
    fid: int = fid_field
    nightdate: date = Field()  # FIXME : why different to obsdate?
    fieldid: int = fieldid_field
    itid: int = Field(ge=0)
    puid: int = Field(ge=0)

    timeutc: datetime = Field()
    AExpTime: float = Field(ge=0)
    expMJD: float = Field(ge=59000)
    airmass: float = Field(ge=1.0)
    shutopen: datetime = Field()
    shutclsd: datetime = Field()
    tempture: float = default_unknown_field
    windspd: float = default_unknown_field
    Dewpoint: float = default_unknown_field
    Humidity: float = default_unknown_field
    Pressure: float = default_unknown_field
    Moonra: float = Field(ge=0.0, le=360.0, default=None)
    Moondec: float = Field(title="Dec (degrees)", ge=-90.0, le=90, default=None)
    Moonillf: float = default_unknown_field
    Moonphas: float = default_unknown_field
    Moonaz: float = default_unknown_field
    Moonalt: float = default_unknown_field
    Sunaz: float = default_unknown_field
    Sunalt: float = default_unknown_field
    Detfirm: str = Field(default="unknown")
    Detsoft: str = Field(default="unknown")
    ra: float = ra_field
    dec: float = dec_field
    altitude: float = alt_field
    azimuth: float = az_field

    def insert_entry(self, returning_key_names=None) -> tuple:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: None
        """
        night = Night(nightdate=self.nightdate)
        if not night.exists():
            night.insert_entry()

        logger.debug(f"puid: {self.puid}")
        if not ProgramsTable().exists(values=self.puid, keys="puid"):
            default_puid = ProgramsTable().select_query(
                compare_values=default_program.progname,
                compare_keys="progname",
                select_keys="puid",
            )[0][0]
            logger.debug(f"Found progid {default_puid}")
            self.puid = default_puid

        return self._insert_entry()

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(values=self.expid, keys="expid")
