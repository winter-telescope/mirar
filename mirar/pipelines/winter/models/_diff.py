"""
Models for the 'diff' table
"""
import os
from typing import ClassVar

from pydantic import Field, field_validator
from sqlalchemy import VARCHAR, Column, ForeignKey, Integer, Sequence
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB
from mirar.pipelines.winter.models._exposures import Exposure
from mirar.pipelines.winter.models.base_model import WinterBase


class DiffsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Diffs table in database
    """

    __tablename__ = "diffs"
    __table_args__ = {"extend_existing": True}

    diffid = Column(
        Integer,
        Sequence(start=1, name="diff_udiffid_seq"),
        autoincrement=True,
        unique=True,
        primary_key=True,
    )

    uexpid: Mapped[int] = mapped_column(ForeignKey("exposures.uexpid"))
    exposure_ids: Mapped["ExposuresTable"] = relationship(back_populates="diff")

    rawid: Mapped[int] = mapped_column(ForeignKey("raws.rawid"))
    raw_ids: Mapped["RawsTable"] = relationship(back_populates="diff")

    stackid: Mapped[int] = mapped_column(ForeignKey("stacks.stackid"))
    stack_id: Mapped["StacksTable"] = relationship(back_populates="diff")

    savepath = Column(VARCHAR(255), unique=True)


class Diff(BaseDB):
    """
    A pydantic model for a diff database entry
    """

    sql_model: ClassVar = DiffsTable

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

        :param savepath: field value
        :return: field value
        """
        assert os.path.exists(savepath)
        return savepath

    @field_validator("uexpid")
    @classmethod
    def validate_expid(cls, uexpid: int) -> int:
        """
        Ensure that expid exists in exposures table

        :param uexpid: field value
        :return: field value
        """
        assert Exposure._exists(keys="uexpid", values=uexpid)
        return uexpid
