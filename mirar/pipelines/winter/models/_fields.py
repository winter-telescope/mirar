"""
Models for the 'fields' table
"""
import time
from typing import ClassVar

from pydantic import Field
from sqlalchemy import REAL, Column, Insert, Integer, Select
from sqlalchemy.orm import Mapped, relationship
from tqdm import tqdm
from wintertoo.data import winter_fields

from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.sqldatabase.base_model import BaseDB, _exists, dec_field, ra_field
from mirar.utils.sql import get_engine

DEFAULT_FIELD = 999999999


class FieldsTable(WinterBase):  # pylint: disable=too-few-public-methods
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


class Fields(BaseDB):
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
    Downloads a field grid (text file) and imports it in chunks into the database

    :return: None
    """

    engine = get_engine(db_name=FieldsTable.db_name)
    if not _exists(Select(FieldsTable), engine=engine):
        chunk = 10000

        winter_fields["fieldid"] = winter_fields["ID"]
        winter_fields["ra"] = winter_fields["RA"]
        winter_fields["dec"] = winter_fields["Dec"]
        winter_fields["ebv"] = winter_fields["Ebv"]
        winter_fields["gall"] = winter_fields["Gal_Long"]
        winter_fields["galb"] = winter_fields["Gal_Lat"]

        keys = list(Fields.__fields__)

        idx = list(range(0, len(winter_fields), chunk)) + [len(winter_fields)]

        for k, i in tqdm(enumerate(idx[:-1]), total=len(idx) - 1):
            j = idx[k + 1]

            res = winter_fields[i:j]

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
