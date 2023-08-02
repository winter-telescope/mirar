"""
Models for the 'candidates' table
"""
import logging
from datetime import date, datetime
from typing import ClassVar

from pydantic import Field
from sqlalchemy import (  # event,
    VARCHAR,
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    Sequence,
)
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB, alt_field, az_field, dec_field, ra_field
from mirar.pipelines.winter.models._fields import FieldsTable, fieldid_field
from mirar.pipelines.winter.models._filters import FiltersTable, fid_field
from mirar.pipelines.winter.models._img_type import ImgTypesTable
from mirar.pipelines.winter.models._nights import Night, NightsTable
from mirar.pipelines.winter.models._programs import ProgramsTable, default_program
from mirar.pipelines.winter.models.base_model import WinterBase

logger = logging.getLogger(__name__)


class CandidatesTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "candidates"
    __table_args__ = {"extend_existing": True}

    candid = Column(
        Integer,
        Sequence(name="candidates_uexpid_seq", start=1, increment=1),
        unique=True,
        autoincrement=True,
    )

    jd = Column(Float, nullable=False)
    diffmaglim = Column(Float, nullable=False)
    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))

    progname: Mapped[str] = mapped_column(ForeignKey("programs.progname"))

    isdiffpos = Column(Boolean, nullable=False)

    fieldid: Mapped[int] = mapped_column(ForeignKey("fields.fieldid"))

    xpos = Column(Float, nullable=True)
    ypos = Column(Float, nullable=True)

    ra = Column(Float)
    dec = Column(Float)
    ra_column_name = "ra"
    dec_column_name = "dec"

    magpsf = Column(Float, nullable=False)
    sigmapsf = Column(Float, nullable=False)
    magap = Column(Float, nullable=False)
    sigmagap = Column(Float, nullable=False)
    magapbig = Column(Float, nullable=False)
    sigmagapbig = Column(Float, nullable=False)

    distnr = Column(Float, nullable=False)
    magnr = Column(Float, nullable=False)
    sigmagnr = Column(Float, nullable=False)
    chinr = Column(Float, nullable=False)
    sharpnr = Column(Float, nullable=False)

    sky = Column(Float, nullable=False)
    magdiff = Column(Float, nullable=False)
    fwhm = Column(Float, nullable=False)
    classtar = Column(Float, nullable=False)
    mindtoedge = Column(Float, nullable=False)
    magfromlim = Column(Float, nullable=False)
    seeratio = Column(Float, nullable=False)
    aimage = Column(Float, nullable=False)
    bimage = Column(Float, nullable=False)
    aimagerat = Column(Float, nullable=False)
    bimagerat = Column(Float, nullable=False)
    elong = Column(Float, nullable=False)
    nneg = Column(Integer, nullable=False)
    nbad = Column(Integer, nullable=False)

    rb = Column(Float, nullable=True)
    rbversion = Column(Float, nullable=True)

    ssdistnr = Column(Float, nullable=True)
    ssmagnr = Column(Float, nullable=True)
    ssnamenr = Column(VARCHAR(40), nullable=True)

    psra1 = Column(Float, nullable=True)
    psdec1 = Column(Float, nullable=True)
    psobjid1 = Column(VARCHAR(40), nullable=True)

    jdstarthist = Column(Float, nullable=False)
    jdendhist = Column(Float, nullable=False)

    scorr = Column(Float, nullable=False)

    tooflag = Column(Boolean, nullable=False)

    nightdate: Mapped[int] = mapped_column(ForeignKey("nights.nightdate"))

    utctime = Column(DateTime(timezone=True))

    ExpTime = Column(Float, nullable=False)
    expMJD = Column(Float, nullable=False)
    airmass = Column(Float)
    tempture = Column(Float, default=-999)
    windspd = Column(Float, default=-999)
    Dewpoint = Column(Float, default=-999)
    Humidity = Column(Float, default=-999)
    Pressure = Column(Float, default=-999)

    Moonaz = Column(Float, default=-999)
    Moonalt = Column(Float, default=-999)
    Sunalt = Column(Float, default=-999)

    altitude = Column(Float)
    azimuth = Column(Float)

    raw: Mapped["RawTable"] = relationship(back_populates="candidate_ids")


default_unknown_field = Field(default=-999)


class Candidate(BaseDB):
    """
    A pydantic model for a candidate database entry
    """

    sql_model: ClassVar = CandidatesTable

    expid: int = Field(ge=0)
    fid: int = fid_field
    nightdate: date = Field()  # FIXME : why different to obsdate?
    fieldid: int = fieldid_field
    itid: int = Field(ge=0)
    progname: str = Field(min_length=1)

    utctime: datetime = Field()
    ExpTime: float = Field(ge=0)
    expMJD: float = Field(ge=59000)

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
        logger.debug(f"Searched for night {self.nightdate}")
        if not night.exists():
            night.insert_entry()

        if not ProgramsTable().exists(values=self.progname, keys="progname"):
            default_progname = ProgramsTable().select_query(
                select_keys="progname",
                compare_values=[default_program.progname],
                compare_keys=["progname"],
            )[0][0]
            self.progname = default_progname

        return self._insert_entry()

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(values=self.expid, keys="expid")
