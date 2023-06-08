"""
Models for the 'field' table
"""
import time
from typing import ClassVar

from pydantic import Field
from sqlalchemy import REAL, Column, Insert, Integer, Select
from sqlalchemy.orm import Mapped, relationship
from tqdm import tqdm
from wintertoo.data import summer_fields

from mirar.pipelines.summer.models.base_model import SummerBase
from mirar.processors.sqldatabase.base_model import BaseDB, _exists, dec_field, ra_field
from mirar.utils.sql import get_engine

DEFAULT_FIELD = 999999999


class FieldsTable(SummerBase):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "fields"

    fieldid = Column(Integer, primary_key=True)
    ra = Column(REAL, nullable=True)
    dec = Column(REAL, nullable=True)
    ebv = Column(REAL, nullable=True)
    gall = Column(REAL, nullable=True)
    galb = Column(REAL, nullable=True)
    exposures: Mapped["ExposuresTable"] = relationship(back_populates="field")


fieldid_field: int = Field(ge=0, default=DEFAULT_FIELD, exclude={-99})


class FieldEntry(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = FieldsTable
    fieldid: int = fieldid_field
    ra: float = ra_field
    dec: float = dec_field
    ebv: float = Field(ge=1)
    gall: float = ra_field
    galb: float = dec_field
    # Need to fix formatting of fields file in wintertoo before including these.


def populate_fields():
    """
    Populates the fields table in the database

    :return: None
    """

    engine = get_engine(db_name=FieldsTable.db_name)
    if not _exists(Select(FieldsTable), engine=engine):
        chunk = 10000

        summer_fields["fieldid"] = summer_fields["ID"]
        summer_fields["ra"] = summer_fields["RA"]
        summer_fields["dec"] = summer_fields["Dec"]
        summer_fields["ebv"] = summer_fields["Ebv"]
        summer_fields["gall"] = summer_fields["Gal_Long"]
        summer_fields["galb"] = summer_fields["Gal_Lat"]

        keys = list(FieldEntry.__fields__)

        idx = list(range(0, len(summer_fields), chunk)) + [len(summer_fields)]

        for k, i in tqdm(enumerate(idx[:-1]), total=len(idx) - 1):
            j = idx[k + 1]

            res = summer_fields[i:j]

            stmt = Insert(FieldsTable).values(
                res[keys].to_dict(orient="records"),
            )

            with engine.connect() as conn:
                with conn.begin():
                    conn.execute(stmt)

            time.sleep(1)

        stmt = Insert(FieldsTable).values(fieldid=DEFAULT_FIELD, ra=None, dec=None)

        with engine.connect() as conn:
            with conn.begin():
                conn.execute(stmt)
