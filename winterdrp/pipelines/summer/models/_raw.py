"""
Models for the 'raw' table
"""
from datetime import date, datetime
from typing import ClassVar

from pydantic import Field
from sqlalchemy import (
    DATE,
    VARCHAR,
    Column,
    DateTime,
    Double,
    Float,
    ForeignKey,
    Integer,
    Sequence,
    event,
)
from sqlalchemy.orm import Mapped, mapped_column, relationship

from winterdrp.pipelines.summer.models._fields import FieldsTable, fieldid_field
from winterdrp.pipelines.summer.models._filters import FiltersTable, fid_field
from winterdrp.pipelines.summer.models._itid import ITIDsTable
from winterdrp.pipelines.summer.models._nights import Nights, NightsTable
from winterdrp.pipelines.summer.models._programs import (
    ProgramsTable,
    default_program,
    program_id_field,
)
from winterdrp.pipelines.summer.models.basemodel import (
    SummerBase,
    alt_field,
    az_field,
    date_field,
    dec_field,
    ra_field,
)
from winterdrp.processors.sqldatabase.basemodel import BaseDB
from winterdrp.utils.sql import create_q3c_extension

db_name = "summertest"


class RawTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "raw"
    __table_args__ = {"extend_existing": True}

    rawid = Column(Double, primary_key=True, autoincrement=False)
    rawcount = Column(Integer, Sequence("raw_rawcount_seq", start=1))
    filtername = Column(VARCHAR(20), unique=True)

    nightid: Mapped[int] = mapped_column(ForeignKey("nights.nightid"))
    night: Mapped["NightsTable"] = relationship(back_populates="raw")

    fieldid: Mapped[int] = mapped_column(ForeignKey("fields.fieldid"))
    field: Mapped["FieldsTable"] = relationship(back_populates="raw")

    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))
    filt: Mapped["FiltersTable"] = relationship(back_populates="raw")

    itid: Mapped[int] = mapped_column(ForeignKey("itid.itid"))
    img_type: Mapped["ITIDsTable"] = relationship(back_populates="raw")

    progID: Mapped[int] = mapped_column(ForeignKey("programs.id"))
    programid: Mapped["ProgramsTable"] = relationship(back_populates="raw")

    expid = Column(Double, unique=True)
    savepath = Column(VARCHAR(255), unique=True)
    obsdate = Column(DATE)
    timeutc = Column(DateTime(timezone=True))
    obsID = Column(Integer)

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
    procflag = Column(Integer, default=0)

    ra_column_name = "ra"
    dec_column_name = "dec"


# @event.listens_for(target=RawTable.__table__, identifier="after_create")
# def raw_q3c(tbl, conn, *args, **kw):
#     create_q3c_extension(
#         conn=conn,
#         __tablename__=RawTable.__tablename__,
#         ra_column_name=RawTable.ra_column_name,
#         dec_column_name=RawTable.dec_column_name,
#     )


default_unknown_field = Field(default=-999)


class Raw(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = RawTable

    rawid: int = Field(ge=0)
    nightid: date = Field()  # FIXME : why different to obsdate?
    fieldid: int = fieldid_field
    fid: int = fid_field

    expid: int = Field()
    savepath: str = Field(min_length=1)
    obsdate: date = date_field
    timeutc: datetime = Field()
    obsID: int = Field()
    itid: int = Field()
    AExpTime: float = Field(ge=0)
    expMJD: float = Field(ge=59000)
    airmass: float = Field(ge=1.0)
    shutopen: datetime = Field()
    shutclsd: datetime = Field()
    altitude: float = alt_field
    azimuth: float = az_field
    ra: float = ra_field
    dec: float = dec_field

    progID: int = program_id_field

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

    # @validator("savepath")
    # @classmethod
    # def validate_time_allocation(cls, field_value: Path):
    #     """
    #     Ensure that path exists
    #
    #     :param field_value: field value
    #     :return: field value
    #     """
    #     assert field_value.exists()
    #     return field_value

    def insert_entry(self, returning_keys=None) -> tuple:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: None
        """
        night = Nights(nightid=self.nightid)
        if not night.exists():
            night.insert_entry()

        if not ProgramsTable().exists(value=self.progID, key="progid"):
            self.progID = default_program.progid

        return self._insert_entry()

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(value=self.rawid, key="rawid")
