"""
Models for the 'raw' table
"""
import os
from typing import ClassVar

from pydantic import Field, validator
from sqlalchemy import VARCHAR, Column, ForeignKey, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship

from winterdrp.pipelines.summer.models._exposures import Exposures
from winterdrp.pipelines.summer.models.basemodel import SummerBase
from winterdrp.processors.sqldatabase.basemodel import BaseDB


class RawTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "raw"
    __table_args__ = {"extend_existing": True}

    rawid = Column(Integer, primary_key=True)

    expid: Mapped[int] = mapped_column(ForeignKey("exposures.expid"))
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

    expid: int = Field(ge=0)
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

    @validator("expid")
    @classmethod
    def validate_expid(cls, field_value: int):
        """
        Ensure that expid exists in exposures table
        Args:
            field_value: expid

        Returns:

        """
        assert Exposures.sql_model().exists(keys="expid", values=field_value)
        return field_value
