"""
Models for the 'diff' table
"""

# pylint: disable=duplicate-code
import os
from typing import ClassVar

from pydantic import Field, field_validator
from sqlalchemy import VARCHAR, Column, Double, ForeignKey, Integer, Sequence
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB
from mirar.pipelines.summer.models._exposures import Exposure
from mirar.pipelines.summer.models.base_model import SummerBase


class DiffTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Diff table in database
    """

    __tablename__ = "diff"
    __table_args__ = {"extend_existing": True}

    udiffid = Column(
        Integer,
        Sequence(start=1, name="diff_udiffid_seq"),
        autoincrement=True,
        unique=True,
    )
    diffid = Column(Double, primary_key=True, autoincrement=False)
    proc: Mapped["ProcTable"] = relationship(back_populates="diff_ids")

    uexpid: Mapped[int] = mapped_column(ForeignKey("exposures.uexpid"))
    exposure_ids: Mapped["ExposuresTable"] = relationship(back_populates="diff")

    savepath = Column(VARCHAR(255), unique=True)


class Diff(BaseDB):
    """
    A pydantic model for a diff database entry
    """

    sql_model: ClassVar = DiffTable

    diffid: int = Field(ge=0)
    uexpid: int = Field(ge=0)
    qid: int = Field(ge=0)
    savepath: str = Field(min_length=1)
    procstatus: int = Field(ge=0, default=0)

    @field_validator("savepath")
    @classmethod
    def validate_savepath(cls, savepath: str) -> str:
        """
        Ensure that path exists

        :param savepath: path
        :return: savepath
        """
        assert os.path.exists(savepath)
        return savepath

    @field_validator("uexpid")
    @classmethod
    def validate_expid(cls, uexpid: int):
        """
        Ensure that expid exists in exposures table

        :param uexpid: field value
        :return: field value
        """
        assert Exposure._exists(keys="uexpid", values=uexpid)
        return uexpid
