"""
Models for the 'raw' table
"""
import os
from typing import ClassVar

from pydantic import Field, validator
from sqlalchemy import VARCHAR, Column, Double, Float, ForeignKey, Integer, Sequence

# event,
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.pipelines.winter.models._exposures import Exposures
from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.sqldatabase.base_model import BaseDB


class RawTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Raw table in database
    """

    __tablename__ = "raw"
    __table_args__ = {"extend_existing": True}

    urawid = Column(
        Integer, Sequence(start=1, name="raw_urawid_seq"), autoincrement=True
    )
    rawid = Column(Double, primary_key=True, unique=True, autoincrement=False)

    uexpid: Mapped[int] = mapped_column(ForeignKey("exposures.uexpid"))
    exposure_ids: Mapped["ExposuresTable"] = relationship(back_populates="raw")

    subdetid: Mapped[int] = mapped_column(ForeignKey("subdets.subdetid"))
    subdets: Mapped["SubdetsTable"] = relationship(back_populates="raw")
    t_roic = Column(Float)

    savepath = Column(VARCHAR(255), unique=True)

    # procstatus = Column(Integer, default=0)
    ustackid: Mapped[int] = mapped_column(ForeignKey("stacks.ustackid"), nullable=True)
    stacks: Mapped["StacksTable"] = relationship(back_populates="raw")

    # proc: Mapped["ProcTable"] = relationship(back_populates="raw_ids")
    astrometry: Mapped["AstrometryStatsTable"] = relationship(
        back_populates="astrom_raw_ids"
    )


class Raw(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = RawTable

    rawid: int = Field(ge=0)
    uexpid: int = Field(ge=0)
    subdetid: int = Field(ge=0)
    savepath: str = Field(min_length=1)
    t_roic: float = Field()
    ustackid: int = Field(ge=0, default=None)

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
        assert Exposures.sql_model().exists(keys="uexpid", values=field_value)
        return field_value
