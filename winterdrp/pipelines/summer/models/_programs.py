"""
Models for the 'program' table
"""
from datetime import date
from typing import ClassVar

from pydantic import BaseModel, Field, validator
from sqlalchemy import CHAR, DATE, REAL, VARCHAR, Column, Integer, Select

from winterdrp.pipelines.summer.models.basemodel import (
    Base,
    BaseDB,
    _exists,
    date_field,
)
from winterdrp.utils.security import generate_key

_LEN_PROG_KEY = 20


class ProgramsTable(Base):  # pylint: disable=too-few-public-methods
    """
    Program table in database
    """

    __tablename__ = "program"

    id = Column(Integer, primary_key=True)
    progname = Column(CHAR(8), unique=True)
    prog_key = Column(CHAR(_LEN_PROG_KEY))
    progid = Column(Integer)
    progtitle = Column(VARCHAR(20))
    piname = Column(VARCHAR(20))
    startdate = Column(DATE)
    enddate = Column(DATE)
    hours_allocated = Column(REAL)
    hours_remaining = Column(REAL)
    basepriority = Column(REAL)


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
    progid: int = Field()
    progtitle: str = Field(min_length=1)
    piname: str = Field(min_length=1)
    startdate: date = date_field
    enddate: date = date_field
    hours_allocated: float = Field(ge=0.0)
    hours_remaining: float = Field(ge=0.0)
    basepriority: float = Field(ge=0.0, example=100.0)

    @validator("enddate")
    def check_date(cls, field_value, values):
        """
        Ensure dates are correctly formatted

        :param value: value
        :return: value
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
        return _exists(
            Select(self.sql_model).where(self.sql_model.progname == self.progname)
        )


default_program = Programs(
    progname="2001A000",
    prog_key=generate_key(_LEN_PROG_KEY),
    piname="HAL",
    progid=1,
    progtitle="Auto-pilot",
    startdate=date(2001, 1, 1),
    enddate=date(3001, 1, 1),
    hours_allocated=0,
    hours_remaining=0,
    basepriority=0,
)
