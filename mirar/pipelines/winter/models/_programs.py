"""
Models for the 'program' table
"""

from datetime import date
from typing import ClassVar

from pydantic import BaseModel, Field, model_validator
from sqlalchemy import CHAR, DATE, REAL, VARCHAR, Column, Integer
from sqlalchemy.orm import Mapped, relationship

from mirar.database.base_model import BaseDB, date_field
from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.utils.security import generate_key

LEN_PROG_KEY = 60
DEFAULT_MAX_PRIORITY = 100.0


class ProgramsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Program table in database
    """

    __tablename__ = "programs"

    puid = Column(Integer, primary_key=True)  # Serial key for counting programs
    progid = Column(Integer)  # identifier for program? 0/1/2/3
    progname = Column(CHAR(8), unique=True)  # e.g. 2022A000
    prog_key = Column(CHAR(LEN_PROG_KEY), unique=True)  # unique hash key for program
    pi_name = Column(VARCHAR(20))  # PI Name
    pi_email = Column(VARCHAR(70))  # PI email
    startdate = Column(DATE)  # Start time of program
    enddate = Column(DATE)  # End time of program
    hours_allocated = Column(REAL)  # Total hours allocated
    hours_used = Column(REAL)  # Total hours used
    maxpriority = Column(REAL, default=DEFAULT_MAX_PRIORITY)  # Base priority

    progtitle = Column(VARCHAR(20), nullable=True)  # Optional 20 char descr. of title

    exposures: Mapped["ExposuresTable"] = relationship(back_populates="program_name")


prog_field: str = Field(min_length=8, max_length=8, example="2020A000")


class ProgramCredentials(BaseModel):
    """
    Program credentials to access a program
    """

    progname: str = prog_field
    prog_key: str = Field(
        min_length=LEN_PROG_KEY,
        max_length=LEN_PROG_KEY,
        description="The auto-generated program hash key",
    )


class Program(BaseDB, ProgramCredentials):
    """
    A pydantic model for a program database entry
    """

    sql_model: ClassVar = ProgramsTable
    progid: int = Field(default=1)
    pi_name: str = Field(min_length=1, example="Hubble")
    pi_email: str = Field(min_length=1, example="someone@institute.com")
    startdate: date = date_field
    enddate: date = date_field
    hours_allocated: float = Field(ge=0.0)
    hours_used: float = Field(ge=0.0)
    maxpriority: float = Field(
        ge=DEFAULT_MAX_PRIORITY,
        default=DEFAULT_MAX_PRIORITY,
        description="Max priority",
    )
    progtitle: str = Field(min_length=1, example="A program title")

    @model_validator(mode="after")
    def check_date(self):
        """
        Ensure dates are correctly formatted

        :return: self
        """
        startdate = self.startdate
        enddate = self.enddate
        assert enddate > startdate
        return self

    @model_validator(mode="after")
    def validate_time_allocation(self):
        """
        Ensure that time remaining has a sensible value

        :return: self
        """
        total_time = self.hours_allocated
        hours_used = self.hours_used
        assert not hours_used > total_time
        assert not hours_used < 0
        return self

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(values=self.progname, keys="progname")


default_program = Program(
    progid=1,
    progname="2001A000",
    prog_key=generate_key(LEN_PROG_KEY),
    pi_name="HAL",
    pi_email="winter@miaow.com",
    progtitle="Auto-pilot",
    startdate=date(2001, 1, 1),
    enddate=date(3001, 1, 1),
    hours_allocated=0,
    hours_used=0,
    maxpriority=DEFAULT_MAX_PRIORITY,
)


def populate_programs():
    """
    Populate the programs table with the default program

    :return: None
    """

    if not default_program.exists():
        default_program.insert_entry()
