"""
Models for the 'candidates' table
"""

import logging
from typing import ClassVar

import pandas as pd
from pydantic import Field
from sqlalchemy import (
    VARCHAR,
    BigInteger,
    Boolean,
    Column,
    Float,
    ForeignKey,
    Integer,
    Sequence,
)
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions import select_from_table
from mirar.pipelines.winter.models._fields import fieldid_field
from mirar.pipelines.winter.models._filters import fid_field
from mirar.pipelines.winter.models._programs import Program, default_program
from mirar.pipelines.winter.models.base_model import WinterBase

logger = logging.getLogger(__name__)

CANDIDATE_PREFIX = "WNTR"
NAME_START = "aaaaa"

MIN_NAME_LENGTH = len(CANDIDATE_PREFIX) + len(NAME_START) + 2


class CandidatesTable(WinterBase):  # pylint: disable=too-few-public-methods
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
    deprecated = Column(Boolean, nullable=False, default=False)
    jd = Column(Float, nullable=False)

    # Image properties

    diffid: Mapped[int] = mapped_column(ForeignKey("diffs.diffid"))  # FIXME
    diff_id: Mapped["DiffsTable"] = relationship(back_populates="candidates")

    stackid: Mapped[int] = mapped_column(ForeignKey("stacks.stackid"))  # FIXME
    stack_id: Mapped["StacksTable"] = relationship(back_populates="candidates")

    fid: Mapped[int] = mapped_column(ForeignKey("filters.fid"))
    exptime = Column(Float, nullable=False)

    progname: Mapped[str] = mapped_column(ForeignKey("programs.progname"))

    isdiffpos = Column(Boolean, nullable=False)

    fieldid: Mapped[int] = mapped_column(ForeignKey("fields.fieldid"))

    # Positional properties

    ra = Column(Float)
    dec = Column(Float)
    ra_column_name = "ra"
    dec_column_name = "dec"

    # Zero point properties

    magzpsci = Column(Float, nullable=True)
    magzpsciunc = Column(Float, nullable=True)

    # Photometry properties

    diffmaglim = Column(Float, nullable=True)

    magpsf = Column(Float, nullable=False)
    sigmapsf = Column(Float, nullable=False)

    chipsf = Column(Float, nullable=True)

    magap = Column(Float, nullable=False)
    sigmagap = Column(Float, nullable=False)

    magapbig = Column(Float, nullable=True)
    sigmagapbig = Column(Float, nullable=True)

    magdiff = Column(Float, nullable=True)
    magfromlim = Column(Float, nullable=True)

    # Diagnostic properties

    distnr = Column(Float, nullable=True)
    magnr = Column(Float, nullable=True)
    sigmagnr = Column(Float, nullable=True)

    xpos = Column(Float, nullable=True)
    ypos = Column(Float, nullable=True)

    sky = Column(Float, nullable=True)
    fwhm = Column(Float, nullable=False)
    mindtoedge = Column(Float, nullable=True)
    seeratio = Column(Float, nullable=True)
    aimage = Column(Float, nullable=False)
    bimage = Column(Float, nullable=False)
    aimagerat = Column(Float, nullable=False)
    bimagerat = Column(Float, nullable=False)
    elong = Column(Float, nullable=False)
    nneg = Column(Integer, nullable=True)
    nbad = Column(Integer, nullable=True)
    sumrat = Column(Float, nullable=True)

    # Diff/ZOGY properties

    dsnrms = Column(Float, nullable=True)
    ssnrms = Column(Float, nullable=True)
    dsdiff = Column(Float, nullable=True)

    scorr = Column(Float, nullable=False)

    # Real/bogus properties

    rb = Column(Float, nullable=True)
    rbversion = Column(Float, nullable=True)

    # Solar system properties

    ssdistnr = Column(Float, nullable=True)
    ssmagnr = Column(Float, nullable=True)
    ssnamenr = Column(VARCHAR(40), nullable=True)

    # ToO

    tooflag = Column(Boolean, nullable=False)

    # Cross-match properties

    # Ps1 properties

    psobjectid1 = Column(BigInteger, nullable=True)
    sgmag1 = Column(Float, nullable=True)
    srmag1 = Column(Float, nullable=True)
    simag1 = Column(Float, nullable=True)
    szmag1 = Column(Float, nullable=True)
    distpsnr1 = Column(Float, nullable=True)
    sgscore1 = Column(Float, nullable=True)

    psobjectid2 = Column(BigInteger, nullable=True)
    sgmag2 = Column(Float, nullable=True)
    srmag2 = Column(Float, nullable=True)
    simag2 = Column(Float, nullable=True)
    szmag2 = Column(Float, nullable=True)
    distpsnr2 = Column(Float, nullable=True)
    sgscore2 = Column(Float, nullable=True)

    psobjectid3 = Column(BigInteger, nullable=True)
    sgmag3 = Column(Float, nullable=True)
    srmag3 = Column(Float, nullable=True)
    simag3 = Column(Float, nullable=True)
    szmag3 = Column(Float, nullable=True)
    distpsnr3 = Column(Float, nullable=True)
    sgscore3 = Column(Float, nullable=True)

    # 2Mass properties

    tmjmag1 = Column(Float, nullable=True)
    tmhmag1 = Column(Float, nullable=True)
    tmkmag1 = Column(Float, nullable=True)
    tmobjectid1 = Column(VARCHAR(40), nullable=True)

    tmjmag2 = Column(Float, nullable=True)
    tmhmag2 = Column(Float, nullable=True)
    tmkmag2 = Column(Float, nullable=True)
    tmobjectid2 = Column(VARCHAR(40), nullable=True)

    tmjmag3 = Column(Float, nullable=True)
    tmhmag3 = Column(Float, nullable=True)
    tmkmag3 = Column(Float, nullable=True)
    tmobjectid3 = Column(VARCHAR(40), nullable=True)

    # Gaia properties

    neargaia = Column(Float, nullable=True)
    neargaiabright = Column(Float, nullable=True)
    maggaia = Column(Float, nullable=True)
    maggaiabright = Column(Float, nullable=True)


class Candidate(BaseDB):
    """
    A pydantic model for a candidate database entry
    """

    sql_model: ClassVar = CandidatesTable

    objectid: str = Field(min_length=MIN_NAME_LENGTH)
    deprecated: bool = Field(default=False)

    jd: float = Field(ge=0)

    diffid: int | None = Field(ge=0, default=None)
    stackid: int = Field(ge=0)

    fid: int = fid_field
    exptime: float = Field(ge=0)

    progname: str = Field(min_length=8, max_length=8, example="2020A000")

    isdiffpos: bool = Field(default=True)

    fieldid: int = fieldid_field

    ra: float = ra_field
    dec: float = dec_field

    magzpsci: float | None = Field(default=None)
    magzpsciunc: float | None = Field(ge=0, default=None)

    diffmaglim: float | None = Field(default=None)

    magpsf: float = Field()
    sigmapsf: float = Field()
    chipsf: float | None = Field(ge=0, default=None)

    magap: float = Field()
    sigmagap: float = Field(ge=0)

    magapbig: float | None = Field(default=None)
    sigmagapbig: float | None = Field(ge=0, default=None)

    magdiff: float | None = Field(default=None)
    magfromlim: float | None = Field(default=None)

    distnr: float | None = Field(ge=0, default=None)
    magnr: float | None = Field(ge=0, default=None)
    sigmagnr: float | None = Field(ge=0, default=None)

    xpos: float | None = Field(ge=0, default=None)
    ypos: float | None = Field(ge=0, default=None)

    sky: float | None = Field(ge=0, default=None)
    fwhm: float = Field(ge=0)
    mindtoedge: float | None = Field(ge=0, default=None)
    seeratio: float | None = Field(ge=0, default=None)
    aimage: float = Field(ge=0)
    bimage: float = Field(ge=0)
    aimagerat: float = Field(ge=0)
    bimagerat: float = Field(ge=0)
    elong: float = Field(ge=0)
    nneg: int | None = Field(ge=0, default=None)
    nbad: int | None = Field(ge=0, default=None)
    sumrat: float | None = Field(default=None)

    dsnrms: float | None = Field(default=None)
    ssnrms: float | None = Field(default=None)
    dsdiff: float | None = Field(default=None)

    scorr: float = Field(ge=0)

    rb: float | None = Field(ge=0, default=None)
    rbversion: float | None = Field(ge=0, default=None)

    ssdistnr: float | None = Field(ge=0, default=None)
    ssmagnr: float | None = Field(ge=0, default=None)
    ssnamenr: str | None = Field(min_length=MIN_NAME_LENGTH, default=None)

    tooflag: bool = Field(default=False)

    psobjectid1: int | None = Field(default=None)
    psobjectid2: int | None = Field(default=None)
    psobjectid3: int | None = Field(default=None)

    sgmag1: float | None = Field(default=None)
    srmag1: float | None = Field(default=None)
    simag1: float | None = Field(default=None)
    szmag1: float | None = Field(default=None)
    sgscore1: float | None = Field(ge=0, default=None)
    distpsnr1: float | None = Field(ge=0, default=None)

    sgmag2: float | None = Field(default=None)
    srmag2: float | None = Field(default=None)
    simag2: float | None = Field(default=None)
    szmag2: float | None = Field(default=None)
    sgscore2: float | None = Field(ge=0, default=None)
    distpsnr2: float | None = Field(ge=0, default=None)

    sgmag3: float | None = Field(default=None)
    srmag3: float | None = Field(default=None)
    simag3: float | None = Field(default=None)
    szmag3: float | None = Field(default=None)
    sgscore3: float | None = Field(ge=0, default=None)
    distpsnr3: float | None = Field(ge=0, default=None)

    tmjmag1: float | None = Field(default=None)
    tmhmag1: float | None = Field(default=None)
    tmkmag1: float | None = Field(default=None)
    tmobjectid1: str | None = Field(default=None)

    tmjmag2: float | None = Field(default=None)
    tmhmag2: float | None = Field(default=None)
    tmkmag2: float | None = Field(default=None)
    tmobjectid2: str | None = Field(default=None)

    tmjmag3: float | None = Field(default=None)
    tmhmag3: float | None = Field(default=None)
    tmkmag3: float | None = Field(default=None)
    tmobjectid3: str | None = Field(default=None)

    neargaia: float | None = Field(ge=0, default=None)
    neargaiabright: float | None = Field(ge=0, default=None)
    maggaia: float | None = Field(default=None)
    maggaiabright: float | None = Field(default=None)

    def insert_entry(self, returning_key_names=None) -> pd.DataFrame:
        """
        Insert the pydantic-ified data into the corresponding sql database

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

        return self._insert_entry(returning_key_names=returning_key_names)
