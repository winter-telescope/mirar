"""
Models for the 'candidates' table
"""
#  pylint: disable=duplicate-code
import logging
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, BigInteger, Boolean, Column, Float, Integer, Sequence
from sqlalchemy.orm import DeclarativeBase

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.database.base_table import BaseTable
from mirar.database.credentials import DB_USER
from mirar.database.setup import setup_database

DB_NAME = "wirc"


logger = logging.getLogger(__name__)

CANDIDATE_PREFIX = "WIRC"
NAME_START = "aaaaa"

MIN_NAME_LENGTH = len(CANDIDATE_PREFIX) + len(NAME_START) + 2


class WircBase(DeclarativeBase, BaseTable):
    """
    Parent class for summer database
    """

    db_name = DB_NAME


class CandidatesTable(WircBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "candidates"
    __table_args__ = {"extend_existing": True}

    # extra avro_path, diff img foreign key etc

    # Core fields
    candid = Column(
        BigInteger,
        Sequence(name="candidates_candid_seq", start=1, increment=1),
        unique=True,
        autoincrement=True,
        primary_key=True,
    )
    objectid = Column(VARCHAR(40), nullable=False, unique=False)

    # Positional properties

    ra = Column(Float)
    dec = Column(Float)
    ra_column_name = "ra"
    dec_column_name = "dec"

    fwhm = Column(Float, nullable=True)

    jd = Column(Float, nullable=False)

    fid = Column(Integer, nullable=False)

    diffimgname = Column(VARCHAR(255), nullable=False)
    sciimgname = Column(VARCHAR(255), nullable=False)
    refimgname = Column(VARCHAR(255), nullable=False)

    magpsf = Column(Float, nullable=True)
    sigmapsf = Column(Float, nullable=True)
    chipsf = Column(Float, nullable=True)

    aimage = Column(Float, nullable=False)
    bimage = Column(Float, nullable=False)
    aimagerat = Column(Float, nullable=False)
    bimagerat = Column(Float, nullable=False)
    elong = Column(Float, nullable=False)

    scorr = Column(Float, nullable=False)

    xpos = Column(Float, nullable=True)
    ypos = Column(Float, nullable=True)

    # Zero point properties

    magzpsci = Column(Float, nullable=True)
    magzpsciunc = Column(Float, nullable=True)

    tmjmag1 = Column(Float, nullable=True)
    tmhmag1 = Column(Float, nullable=True)
    tmkmag1 = Column(Float, nullable=True)
    tmobjectid1 = Column(VARCHAR(25), nullable=True)
    isdiffpos = Column(Boolean, nullable=False)


class Candidate(BaseDB):
    """
    A pydantic model for a candidate database entry
    """

    sql_model: ClassVar = CandidatesTable

    objectid: str = Field(min_length=MIN_NAME_LENGTH)

    ra: float = ra_field
    dec: float = dec_field

    fwhm: float = Field(ge=0)

    jd: float = Field(ge=0)

    fid: int = Field(ge=0)

    diffimgname: str | None = Field(max_length=255, default=None)
    sciimgname: str | None = Field(max_length=255, default=None)
    refimgname: str | None = Field(max_length=255, default=None)

    magpsf: float = Field()
    sigmapsf: float = Field(ge=0)
    chipsf: float | None = Field(ge=0, default=None)

    aimage: float = Field(ge=0)
    bimage: float = Field(ge=0)
    aimagerat: float = Field(ge=0)
    bimagerat: float = Field(ge=0)
    elong: float = Field(ge=0)

    scorr: float = Field(ge=0)

    xpos: float | None = Field(ge=0, default=None)
    ypos: float | None = Field(ge=0, default=None)

    magzpsci: float | None = Field(default=None)
    magzpsciunc: float | None = Field(ge=0, default=None)

    tmjmag1: float | None = Field(default=None)
    tmhmag1: float | None = Field(default=None)
    tmkmag1: float | None = Field(default=None)
    tmobjectid1: str | None = Field(default=None)

    isdiffpos: bool = Field(default=True)


if DB_USER is not None:
    setup_database(WircBase)
