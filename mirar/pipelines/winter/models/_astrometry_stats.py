"""
Models for the 'exposures' table
"""
import logging
from typing import ClassVar

from pydantic import Field
from sqlalchemy import Column, Float, ForeignKey  # event,
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.pipelines.winter.models._raw import RawTable
from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.sqldatabase.base_model import BaseDB, dec_field, ra_field

# from mirar.utils.sql import create_q3c_extension


logger = logging.getLogger(__name__)


class AstrometryStatsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Astrometry stats table in database
    """

    __tablename__ = "astrometry_stats"
    __table_args__ = {"extend_existing": True}

    rawid = mapped_column(ForeignKey("raw.rawid"), primary_key=True)
    astrom_raw_ids: Mapped["RawTable"] = relationship(back_populates="astrometry")
    crval1 = Column(Float)
    crval2 = Column(Float)
    crpix1 = Column(Float)
    crpix2 = Column(Float)
    cd1_1 = Column(Float)
    cd1_2 = Column(Float)
    cd2_1 = Column(Float)
    cd2_2 = Column(Float)
    astunc = Column(Float)
    astfield = Column(Float)
    fwhm_med = Column(Float)
    fwhm_std = Column(Float)
    fwhm_pix = Column(Float)

    ra_column_name = "crval1"
    dec_column_name = "crval2"


# @event.listens_for(target=RawTable.__table__, identifier="after_create")
# def raw_q3c(tbl, conn, *args, **kw):
#     create_q3c_extension(
#         conn=conn,
#         __tablename__=RawTable.__tablename__,
#         ra_column_name=RawTable.ra_column_name,
#         dec_column_name=RawTable.dec_column_name,
#     )


default_unknown_field = Field(default=-999)


class AstrometryStats(BaseDB):
    """
    A pydantic model for a raw database entry
    """

    sql_model: ClassVar = AstrometryStatsTable

    rawid: float = default_unknown_field
    crval1: float = ra_field
    crval2: float = dec_field
    crpix1: float = default_unknown_field
    crpix2: float = default_unknown_field
    cd1_1: float = default_unknown_field
    cd1_2: float = default_unknown_field
    cd2_1: float = default_unknown_field
    cd2_2: float = default_unknown_field
    astunc: float = default_unknown_field
    astfield: float = default_unknown_field
    fwhm_med: float = default_unknown_field
    fwhm_std: float = default_unknown_field
    fwhm_pix: float = default_unknown_field

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(values=self.rawid, keys="rawid")
