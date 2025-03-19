"""
Models for the 'exposures' table
"""

import logging
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Float, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.pipelines.winter.models._raw import RawsTable
from mirar.pipelines.winter.models.base_model import WinterBase

logger = logging.getLogger(__name__)


class FirstPassAstrometryStatsTable(
    WinterBase
):  # pylint: disable=too-few-public-methods
    """
    Astrometry stats table in database
    """

    __tablename__ = "fp_astrometry_stats"
    __table_args__ = {"extend_existing": True}

    rawid = mapped_column(ForeignKey("raws.rawid"), primary_key=True, unique=True)
    fp_astrom_raw_ids: Mapped["RawsTable"] = relationship(
        back_populates="fp_astrometry"
    )
    savepath = Column(VARCHAR(255), unique=True)
    crval1 = Column(Float)
    crval2 = Column(Float)
    crpix1 = Column(Float)
    crpix2 = Column(Float)
    cd1_1 = Column(Float)
    cd1_2 = Column(Float)
    cd2_1 = Column(Float)
    cd2_2 = Column(Float)
    astunc = Column(Float)
    astunc95 = Column(Float)
    astfield = Column(Float)
    fwhm_med = Column(Float)
    fwhm_std = Column(Float)
    fwhm_pix = Column(Float)
    astirms1 = Column(Float)
    astirms2 = Column(Float)
    astrrms1 = Column(Float)
    astrrms2 = Column(Float)

    ra_column_name = "crval1"
    dec_column_name = "crval2"


default_unknown_field = Field(default=-999)


class FirstPassAstrometryStat(BaseDB):
    """
    A pydantic model for an astrometry stats database entry
    """

    sql_model: ClassVar = FirstPassAstrometryStatsTable

    rawid: int = Field(ge=0)
    savepath: str = Column(VARCHAR(255), unique=True)
    crval1: float = ra_field
    crval2: float = dec_field
    crpix1: float = default_unknown_field
    crpix2: float = default_unknown_field
    cd1_1: float = default_unknown_field
    cd1_2: float = default_unknown_field
    cd2_1: float = default_unknown_field
    cd2_2: float = default_unknown_field
    astunc: float = default_unknown_field
    astunc95: float = default_unknown_field
    astfield: float = default_unknown_field
    fwhm_med: float = default_unknown_field
    fwhm_std: float = default_unknown_field
    fwhm_pix: float = default_unknown_field
    astirms1: float = default_unknown_field
    astirms2: float = default_unknown_field
    astrrms1: float = default_unknown_field
    astrrms2: float = default_unknown_field

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(values=self.rawid, keys="rawid")
