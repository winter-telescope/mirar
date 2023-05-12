"""
Models for the 'programs' table
"""
import os
from datetime import date
from typing import ClassVar

import pandas as pd
from pydantic import BaseModel, Field, validator
from sqlalchemy import CHAR, DATE, REAL, VARCHAR, Column, Integer
from sqlalchemy.orm import Mapped, relationship

from winterdrp.pipelines.winter.models.base_model import WinterBase
from winterdrp.processors.sqldatabase.base_model import BaseDB, date_field
from winterdrp.utils.security import generate_key

_LEN_PROG_KEY = 20
program_id_field: int = Field(default=1)


class ProgramsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Program table in database
    """

    __tablename__ = "programs"

    puid = Column(Integer, primary_key=True)  # Serial key for counting programs
    progid = Column(Integer)  # identifier for program? 0/1/2/3
    progname = Column(CHAR(8), unique=True)  # e.g. 2022A000
    prog_key = Column(CHAR(_LEN_PROG_KEY), unique=True)  # unique hash key for program
    pi_name = Column(VARCHAR(20))  # PI Name
    pi_email = Column(VARCHAR(30))  # PI email
    startdate = Column(DATE)  # Start time of program
    enddate = Column(DATE)  # End time of program
    hours_allocated = Column(REAL)  # Total hours allocated
    hours_remaining = Column(REAL)  # Total hours remaining
    basepriority = Column(REAL, default=0)  # Base priority

    progtitle = Column(VARCHAR(20), nullable=True)  # Optional 20 char descr. of title

    exposures: Mapped["ExposuresTable"] = relationship(back_populates="program_uid")


class ProgramCredentials(BaseModel):
    """
    Program credentials to access a program
    """

    progname: str = Field(min_length=8, max_length=8, example="2020A000")
    prog_key: str = Field(min_length=_LEN_PROG_KEY, max_length=_LEN_PROG_KEY)


class Programs(BaseDB, ProgramCredentials):
    """
    A pydantic model for a program database entry
    """

    sql_model: ClassVar = ProgramsTable
    progid: int = program_id_field
    progname: str = Field(min_length=1)
    prog_key: str = Field(min_length=1)
    pi_name: str = Field(min_length=1)
    pi_email: str = Field(min_length=1)
    startdate: date = date_field
    enddate: date = date_field
    hours_allocated: float = Field(ge=0.0)
    hours_remaining: float = Field(ge=0.0)
    basepriority: float = Field(ge=0.0, example=100.0)
    progtitle: str = Field(min_length=1)

    @validator("enddate")
    @classmethod
    def check_date(cls, field_value, values):
        """
        Ensure dates are correctly formatted

        """
        startdate = values["startdate"]
        assert field_value > startdate
        return field_value

    @validator("hours_remaining")
    @classmethod
    def validate_time_allocation(cls, field_value, values):
        """
        Ensure that time remaining has a sensible value

        :param field_value: field value
        :param values: values
        :return: field value
        """
        total_time = values["hours_allocated"]
        assert not field_value > total_time
        assert not field_value < 0.0
        return field_value

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(values=self.progname, keys="progname")


default_program = Programs(
    progid=1,
    progname="2001A000",
    prog_key=generate_key(_LEN_PROG_KEY),
    pi_name="HAL",
    pi_email="winter@miaow.com",
    progtitle="Auto-pilot",
    startdate=date(2001, 1, 1),
    enddate=date(3001, 1, 1),
    hours_allocated=0,
    hours_remaining=0,
    basepriority=0,
)


def populate_programs():
    if not default_program.exists():
        default_program.insert_entry()

    programs_file = os.getenv("WINTER_PROGRAMS_FILE", None)
    if programs_file is not None:
        programs_table = pd.read_csv(programs_file)
        for ind in range(len(programs_table)):
            row = programs_table.iloc[ind]
            prog_key = generate_key(_LEN_PROG_KEY)
            program = Programs(prog_key=prog_key, **row)

            if not program.exists():
                _ = program.insert_entry(returning_key_names="prog_key")
