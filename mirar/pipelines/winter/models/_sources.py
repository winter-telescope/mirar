"""
Models for the 'sources' table
"""

import logging
from datetime import datetime
from typing import ClassVar, List

from pydantic import Field
from sqlalchemy import (
    VARCHAR,
    BigInteger,
    Boolean,
    Column,
    DateTime,
    Float,
    Integer,
    Sequence,
)
from sqlalchemy.orm import Mapped, relationship

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.pipelines.winter.models.base_model import WinterBase

logger = logging.getLogger(__name__)

SOURCE_PREFIX = "WNTR"
NAME_START = "aaaaa"

MIN_NAME_LENGTH = len(SOURCE_PREFIX) + len(NAME_START) + 2


class SourcesTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Sources table in database
    """

    __tablename__ = "sources"
    __table_args__ = {"extend_existing": True}

    # extra avro_path, diff img foreign key etc

    # Core fields
    sourceid = Column(
        BigInteger,
        Sequence(name="sources_candid_seq", start=1, increment=1),
        unique=True,
        autoincrement=True,
        primary_key=True,
    )
    objectid = Column(VARCHAR(40), nullable=False, unique=True)
    deprecated = Column(Boolean, nullable=False, default=False)

    # Positional properties

    average_ra = Column(Float)
    average_dec = Column(Float)
    ra_column_name = "average_ra"
    dec_column_name = "average_dec"

    ndet = Column(Integer, nullable=False, default=1)
    first_det_utc = Column(DateTime(timezone=True))
    latest_det_utc = Column(DateTime(timezone=True))

    candidates: Mapped[List["CandidatesTable"]] = relationship(back_populates="source")


class Source(BaseDB):
    """
    A pydantic model for a source database entry
    """

    sql_model: ClassVar = SourcesTable

    objectid: str = Field(min_length=MIN_NAME_LENGTH)
    deprecated: bool = Field(default=False)

    average_ra: float = ra_field
    average_dec: float = dec_field

    first_det_utc: datetime = Field(description="UTC of first detection")
    latest_det_utc: datetime = Field(description="UTC of latest detection")

    ndet: int = Field(ge=1, description="Number of detections", default=1)
