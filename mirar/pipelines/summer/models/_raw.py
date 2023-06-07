"""
Models for the 'raw' table
"""
import os
from typing import ClassVar

from pydantic import Field, validator
from sqlalchemy import VARCHAR, Column, Double, ForeignKey, Integer, Sequence
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.pipelines.summer.models._exposures import Exposure
from mirar.pipelines.summer.models.base_model import SummerBase
from mirar.processors.sqldatabase.base_model import BaseDB


class RawTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "raw"
    __table_args__ = {"extend_existing": True}

    urawid = Column(
        Integer,
        Sequence(start=1, name="raw_urawid_seq"),
        autoincrement=True,
        unique=True,
    )
    rawid = Column(Double, primary_key=True, autoincrement=False)
    proc: Mapped["ProcTable"] = relationship(back_populates="raw_ids")

    uexpid: Mapped[int] = mapped_column(ForeignKey("exposures.uexpid"))
    exposure_ids: Mapped["ExposuresTable"] = relationship(back_populates="raw")

    qid: Mapped[int] = mapped_column(ForeignKey("subdets.qid"))
    subdets: Mapped["SubdetsTable"] = relationship(back_populates="raw")

    savepath = Column(VARCHAR(255), unique=True)

    procstatus = Column(Integer, default=0)


class Raw(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = RawTable

    rawid: int = Field(ge=0)
    uexpid: int = Field(ge=0)
    qid: int = Field(ge=0)
    savepath: str = Field(min_length=1)
    procstatus: int = Field(ge=0, default=0)

    @validator("savepath")
    @classmethod
    def validate_savepath(cls, field_value: str):
        """
        Ensure that path exists

        :param field_value: field value
        :return: field value
        """
        assert os.path.exists(field_value)
        return field_value

    @validator("uexpid")
    @classmethod
    def validate_expid(cls, field_value: int):
        """
        Ensure that expid exists in exposures table
        Args:
            field_value: expid

        Returns:

        """
        assert Exposure.sql_model().exists(keys="uexpid", values=field_value)
        return field_value
